// Rough throughput benchmark: cnpy_bench [MiB] [directory]
// Numbers include the page cache, i.e. they measure library overhead rather than disk speed.
#include <cnpy/cnpy.hpp>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <numeric>
#include <string>
#include <vector>

namespace {

double best_of(int runs, const std::function<void()>& f) {
    double best = 1e300;
    for (int r = 0; r < runs; ++r) {
        const auto t0 = std::chrono::steady_clock::now();
        f();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double>(t1 - t0).count());
    }
    return best;
}

void report(const char* what, double seconds, double mib) {
    std::printf("%-34s %9.1f ms %9.0f MiB/s\n", what, seconds * 1e3, mib / seconds);
}

} // namespace

int main(int argc, char** argv) {
    const std::size_t mib = argc > 1 ? std::strtoul(argv[1], nullptr, 10) : 256;
    const std::filesystem::path dir = argc > 2 ? argv[2] : std::filesystem::temp_directory_path();
    const std::size_t n = mib * 1024 * 1024 / sizeof(double);
    const auto npy = dir / "cnpy_bench.npy";
    const auto npz = dir / "cnpy_bench.npz";
    const auto npzc = dir / "cnpy_bench_c.npz";

    std::vector<double> data(n);
    for (std::size_t i = 0; i < n; ++i)
        data[i] = static_cast<double>(i % 1000) * 0.5; // compressible
    const auto m = static_cast<double>(mib);
    volatile double sink = 0;

    report("save_npy", best_of(3, [&] { cnpy::save_npy(npy, data, {n}); }), m);
    report("load_npy", best_of(3, [&] { sink = cnpy::load_npy(npy).data<double>()[n / 2]; }), m);
    report("map_npy (open only)", best_of(3, [&] { sink = static_cast<double>(cnpy::map_npy(npy).shape()[0]); }), m);
    report("map_npy + sum",
           best_of(3,
                   [&] {
                       const auto a = cnpy::map_npy(npy);
                       const auto d = a.data<double>();
                       sink = std::accumulate(d.begin(), d.end(), 0.0);
                   }),
           m);
    report("append_npy (2 x half)",
           best_of(3,
                   [&] {
                       std::filesystem::remove(npy);
                       cnpy::append_npy(npy, std::span(data).first(n / 2));
                       cnpy::append_npy(npy, std::span(data).subspan(n / 2));
                   }),
           m);
    report("save_npz (stored)", best_of(3, [&] { cnpy::save_npz(npz, "x", data); }), m);
    report("load_npz (stored, crc)", best_of(3, [&] { sink = cnpy::load_npz(npz, "x").data<double>()[1]; }), m);
    cnpy::npz_read_options no_crc;
    no_crc.verify_crc = false;
    report("load_npz (stored, no crc)", best_of(3, [&] { sink = cnpy::load_npz(npz, "x", no_crc).data<double>()[1]; }),
           m);
    report(
        "save_npz (deflate level 1)",
        best_of(1, [&] { cnpy::save_npz(npzc, "x", data, cnpy::npz_mode::truncate, cnpy::compression::deflate(1)); }),
        m);
    report("load_npz (deflate)", best_of(3, [&] { sink = cnpy::load_npz(npzc, "x").data<double>()[1]; }), m);

    std::filesystem::remove(npy);
    std::filesystem::remove(npz);
    std::filesystem::remove(npzc);
    (void)sink;
    return 0;
}

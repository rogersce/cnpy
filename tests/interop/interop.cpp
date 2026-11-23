// Companion program for interop.py.
//
//   cnpy_interop read  <dir>   verify files written by NumPy
//   cnpy_interop write <dir>   write files for NumPy to verify
//
// Values follow the formula in ../helpers.hpp, mirrored in interop.py.
#include "../helpers.hpp"

#include <complex>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <tuple>

namespace {

template <class T>
struct named {
    using type = T;
    const char* name;
};

const auto types = std::make_tuple(
    named<bool>{"bool"}, named<std::int8_t>{"int8"}, named<std::uint8_t>{"uint8"}, named<std::int16_t>{"int16"},
    named<std::uint16_t>{"uint16"}, named<std::int32_t>{"int32"}, named<std::uint32_t>{"uint32"},
    named<std::int64_t>{"int64"}, named<std::uint64_t>{"uint64"}, named<float>{"float32"}, named<double>{"float64"},
    named<std::complex<float>>{"complex64"}, named<std::complex<double>>{"complex128"});

template <class F>
void for_each_type(F&& f) {
    std::apply([&](auto... t) { (f(t), ...); }, types);
}

int failures = 0;

void fail(const std::string& what) {
    std::cerr << "FAIL: " << what << '\n';
    ++failures;
}

/// Checks a (3, 4, 5) array whose logical C-order element n equals value(n), in either memory order.
template <class T>
void check_3d(const cnpy::ndarray& a, const std::string& what, bool fortran) {
    if (a.shape() != cnpy::shape_t{3, 4, 5} || a.fortran_order() != fortran || !a.dtype().holds<T>())
        return fail(what + ": wrong shape/order/dtype " + a.dtype().str());
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 4; ++j)
            for (std::size_t k = 0; k < 5; ++k)
                if (a.at<T>(i, j, k) != test::value<T>(i * 20 + j * 5 + k)) return fail(what + ": wrong value");
}

template <class T>
void check_flat(const cnpy::ndarray& a, const std::string& what, const cnpy::shape_t& shape) {
    if (a.shape() != shape || !a.dtype().holds<T>()) return fail(what + ": wrong shape/dtype " + a.dtype().str());
    if (!test::same(a.data<T>(), test::values<T>(a.size()))) fail(what + ": wrong values");
}

template <class F>
void guarded(const std::string& what, F&& f) {
    try {
        f();
    } catch (const std::exception& e) {
        fail(what + ": " + e.what());
    }
}

void read_all(const fs::path& dir) {
    for_each_type([&](auto t) {
        using T = typename decltype(t)::type;
        const std::string n = t.name;
        guarded(n + "_c", [&] { check_3d<T>(cnpy::load_npy(dir / (n + "_c.npy")), n + "_c", false); });
        guarded(n + "_f", [&] { check_3d<T>(cnpy::load_npy(dir / (n + "_f.npy")), n + "_f", true); });
        guarded(n + "_be", [&] { check_3d<T>(cnpy::load_npy(dir / (n + "_be.npy")), n + "_be", false); });
        guarded(n + "_mmap", [&] { check_3d<T>(cnpy::map_npy(dir / (n + "_c.npy")), n + "_mmap", false); });
        guarded(n + "_1d", [&] { check_flat<T>(cnpy::load_npy(dir / (n + "_1d.npy")), n + "_1d", {17}); });
        guarded(n + "_0d", [&] { check_flat<T>(cnpy::load_npy(dir / (n + "_0d.npy")), n + "_0d", {}); });
        guarded(n + "_empty", [&] { check_flat<T>(cnpy::load_npy(dir / (n + "_empty.npy")), n + "_empty", {0, 3}); });
        guarded(n + "_v2", [&] { check_flat<T>(cnpy::load_npy(dir / (n + "_v2.npy")), n + "_v2", {2, 3}); });
        guarded(n + "_v3", [&] { check_flat<T>(cnpy::load_npy(dir / (n + "_v3.npy")), n + "_v3", {2, 3}); });
    });

    for (const char* archive : {"savez.npz", "savez_compressed.npz"}) {
        for (const bool mapped : {false, true}) {
            const std::string what = std::string(archive) + (mapped ? " (mmap)" : "");
            guarded(what, [&] {
                cnpy::npz_read_options opt;
                opt.memory_map = mapped;
                const cnpy::npz_reader r(dir / archive, opt);
                std::size_t expected = 0;
                for_each_type([&](auto t) {
                    using T = typename decltype(t)::type;
                    ++expected;
                    check_3d<T>(r.load(t.name), what + "/" + t.name, false);
                });
                if (r.size() != expected) fail(what + ": wrong member count");
            });
        }
    }

    guarded("strings", [&] {
        const auto a = cnpy::load_npy(dir / "strings.npy");
        if (a.dtype().str() != std::string(1, static_cast<char>(cnpy::native_byte_order)) + "U5" ||
            a.shape() != cnpy::shape_t{2} || a.nbytes() != 40)
            fail("strings: wrong dtype/shape " + a.dtype().str());
        const auto b = a.bytes();
        if (std::to_integer<char>(b[0]) != 'h' && std::to_integer<char>(b[3]) != 'h') fail("strings: wrong content");
    });
    guarded("structured", [&] {
        try {
            (void)cnpy::load_npy(dir / "structured.npy");
            fail("structured: expected format_error");
        } catch (const cnpy::format_error&) {
        }
    });
}

void write_all(const fs::path& dir) {
    cnpy::npz_writer stored(dir / "cnpy_savez.npz");
    cnpy::npz_writer packed(dir / "cnpy_savez_compressed.npz");
    for_each_type([&](auto t) {
        using T = typename decltype(t)::type;
        const std::string n = t.name;
        const auto v = test::values<T>(60);
        cnpy::save_npy(dir / (n + "_c.npy"), v, {3, 4, 5});
        cnpy::save_npy(dir / (n + "_1d.npy"), test::values<T>(17));
        cnpy::save_npy(dir / (n + "_0d.npy"), test::values<T>(1), {});
        cnpy::save_npy(dir / (n + "_empty.npy"), test::vec<T>(0), {0, 3});

        // Fortran file: same logical array, column-major memory.
        test::vec<T> f(60);
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 4; ++j)
                for (std::size_t k = 0; k < 5; ++k)
                    f[i + 3 * j + 12 * k] = v[i * 20 + j * 5 + k];
        cnpy::save_npy(dir / (n + "_f.npy"), f, {3, 4, 5}, true);

        // Built by appending 1 + 2 + ... rows: must equal the full (3, 4, 5) array.
        cnpy::append_npy(dir / (n + "_appended.npy"), std::span(v).first(20), {1, 4, 5});
        cnpy::append_npy(dir / (n + "_appended.npy"), std::span(v).subspan(20), {2, 4, 5});

        stored.add(n, v, {3, 4, 5});
        packed.add(n, v, {3, 4, 5}, cnpy::compression::deflate());
    });
    stored.close();
    packed.close();

    // Mixed archive built in three sessions.
    cnpy::save_npz(dir / "cnpy_appended.npz", "a", test::values<double>(10));
    cnpy::save_npz(dir / "cnpy_appended.npz", "b", test::values<std::int32_t>(1000), cnpy::npz_mode::append,
                   cnpy::compression::deflate(9));
    cnpy::save_npz(dir / "cnpy_appended.npz", "c", test::values<float>(6), {2, 3}, cnpy::npz_mode::append);

    // More than 65535 members: requires zip64 end-of-central-directory records.
    cnpy::npz_writer many(dir / "cnpy_many.npz");
    for (std::uint32_t i = 0; i < 65536 + 5; ++i)
        many.add("m" + std::to_string(i), std::span(&i, 1), cnpy::shape_t{});
    many.close();
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: cnpy_interop read|write <dir>\n";
        return 2;
    }
    const std::string mode = argv[1];
    const fs::path dir = argv[2];
    if (mode == "read") {
        read_all(dir);
    } else if (mode == "write") {
        guarded("write", [&] { write_all(dir); });
    } else {
        std::cerr << "unknown mode " << mode << '\n';
        return 2;
    }
    if (failures != 0) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "cnpy_interop " << mode << ": all checks passed\n";
    return 0;
}

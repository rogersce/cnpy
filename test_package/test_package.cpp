#include <cnpy/cnpy.hpp>

#include <cstdio>
#include <vector>

int main() {
    const auto path = std::filesystem::temp_directory_path() / "cnpy_test_package.npz";
    cnpy::save_npz(path, "x", std::vector<double>{1.0, 2.0, 3.0}, cnpy::npz_mode::truncate,
                   cnpy::compression::deflate());
    const double last = cnpy::load_npz(path, "x").data<double>()[2];
    std::filesystem::remove(path);
    std::printf("cnpy test_package: %s\n", last == 3.0 ? "ok" : "FAILED");
    return last == 3.0 ? 0 : 1;
}

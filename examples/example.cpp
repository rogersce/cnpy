// Tour of the cnpy API. Writes a few files into the current directory.
#include <cnpy/cnpy.hpp>

#include <complex>
#include <iostream>
#include <vector>

int main() try {
    constexpr std::size_t nx = 128, ny = 64, nz = 32;

    std::vector<std::complex<double>> data(nx * ny * nz);
    for (std::size_t i = 0; i < data.size(); ++i)
        data[i] = {static_cast<double>(i), -static_cast<double>(i)};

    // .npy: write, read back, access typed data
    cnpy::save_npy("arr1.npy", data, {nz, ny, nx});
    const cnpy::ndarray arr = cnpy::load_npy("arr1.npy");
    const auto values = arr.data<std::complex<double>>(); // std::span, throws on dtype mismatch
    std::cout << "arr1.npy: dtype " << arr.dtype().str() << ", " << arr.ndim() << " dims, " << arr.size()
              << " elements, last = " << values.back() << ", at(1,2,3) = " << arr.at<std::complex<double>>(1, 2, 3)
              << '\n';

    // Grow the file along axis 0 -> shape (2*nz, ny, nx)
    cnpy::append_npy("arr1.npy", data, {nz, ny, nx});

    // Memory-map instead of reading: no copy, pages are loaded on demand
    const cnpy::ndarray mapped = cnpy::map_npy("arr1.npy");
    std::cout << "mapped shape[0] = " << mapped.shape()[0] << '\n';

    // .npz: several arrays in one archive, optionally compressed
    const double scalar = 1.2;
    cnpy::npz_writer npz("out.npz");
    npz.add("scalar", std::span(&scalar, 1), cnpy::shape_t{}); // 0-d array
    npz.add("vec", std::vector<int>{1, 2, 3});
    npz.add("arr1", data, {nz, ny, nx}, cnpy::compression::deflate());
    npz.close();

    const cnpy::npz_reader reader("out.npz");
    for (const auto& name : reader.names())
        std::cout << "out.npz contains '" << name << "'\n";
    std::cout << "scalar = " << reader.load("scalar").at<double>() << '\n';

    const cnpy::npz_t all = cnpy::load_npz("out.npz");
    std::cout << "vec[2] = " << all.at("vec").data<int>()[2] << '\n';
    return 0;
} catch (const cnpy::error& e) {
    std::cerr << "cnpy error: " << e.what() << '\n';
    return 1;
}

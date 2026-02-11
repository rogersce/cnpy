#include <cnpy/cnpy.h>

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

int main() {
    const std::string file = "test_npy_roundtrip.npy";
    std::vector<double> data = {1.25, 2.5, 3.75, 5.0};

    cnpy::npy_save(file, data.data(), {data.size()}, "w");

    cnpy::NpyArray arr = cnpy::npy_load(file);
    const double* loaded = arr.data<double>();

    assert(arr.word_size == sizeof(double));
    assert(arr.shape.size() == 1);
    assert(arr.shape[0] == data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        assert(loaded[i] == data[i]);
    }

    std::remove(file.c_str());
    return 0;
}

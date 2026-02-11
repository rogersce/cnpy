#include <cnpy/cnpy.h>

#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

int main() {
#if CNPY_HAS_ZLIB
    const std::string file = "test_npz_roundtrip.npz";
    std::vector<int> values = {10, 20, 30, 40};
    double scalar = 42.5;

    cnpy::npz_save(file, "values", values.data(), {values.size()}, "w");
    cnpy::npz_save(file, "scalar", &scalar, {1}, "a");

    cnpy::NpyArray loaded_values = cnpy::npz_load(file, "values");
    const int* loaded_data = loaded_values.data<int>();

    assert(loaded_values.word_size == sizeof(int));
    assert(loaded_values.shape.size() == 1);
    assert(loaded_values.shape[0] == values.size());
    for (size_t i = 0; i < values.size(); ++i) {
        assert(loaded_data[i] == values[i]);
    }

    cnpy::npz_t all = cnpy::npz_load(file);
    cnpy::NpyArray loaded_scalar = all["scalar"];
    const double* scalar_ptr = loaded_scalar.data<double>();
    assert(loaded_scalar.shape.size() == 1);
    assert(loaded_scalar.shape[0] == 1);
    assert(*scalar_ptr == scalar);

    std::remove(file.c_str());
#endif
    return 0;
}

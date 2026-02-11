# Purpose

NumPy offers the `save` method for easy saving of arrays into `.npy` and `savez` for zipping multiple `.npy` arrays into a single `.npz`.

`cnpy` lets you read and write these formats directly in C++.

The motivation comes from scientific programming where large amounts of data are generated in C++ and analyzed in Python.

# Project layout

- `include/cnpy/cnpy.h`: public header (header-only library)
- `examples/example1.cpp`: usage example
- `tools/`: helper conversion scripts (`mat2npz`, `npy2mat`, `npz2mat`)

# Installation

Default installation directory is `/usr/local`.
To specify a different directory, add `-DCMAKE_INSTALL_PREFIX=/path/to/install/dir` to the CMake command.

1. Install [CMake](https://www.cmake.org)
2. Create a build directory (for example, `$HOME/build`)
3. `cd $HOME/build`
4. `cmake /path/to/cnpy`
5. `make`
6. `make install`

# Using

Include the library in your code:

```cpp
#include <cnpy/cnpy.h>
```

Compile your program:

```bash
# .npy + .npz (requires zlib)
g++ -o mycode mycode.cpp -I/path/to/install/dir/include -lz --std=c++11

# .npy only
g++ -o mycode mycode.cpp -I/path/to/install/dir/include -DCNPY_HAS_ZLIB=0 --std=c++11
```

If you build via CMake, `.npz` support is controlled by `CNPY_WITH_ZLIB` (default `ON`). Configure with `-DCNPY_WITH_ZLIB=OFF` for `.npy`-only builds.

# API overview

There are two functions for writing data:
- `npy_save`
- `npz_save` *(available when `CNPY_HAS_ZLIB=1`)*

There are three functions for reading:
- `npy_load` loads a `.npy`
- `npz_load(fname)` loads a `.npz` into a map of `NpyArray` *(available when `CNPY_HAS_ZLIB=1`)*
- `npz_load(fname,varname)` loads one variable from a `.npz` *(available when `CNPY_HAS_ZLIB=1`)*

Loaded data is represented by:

```c++
struct NpyArray {
    std::vector<size_t> shape;
    size_t word_size;
    template<typename T> T* data();
};
```

See `examples/example1.cpp` for a complete workflow.

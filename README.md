# cnpy

[![CI](https://github.com/lelapinagile/cnpy/actions/workflows/ci.yml/badge.svg)](https://github.com/lelapinagile/cnpy/actions/workflows/ci.yml)

Read and write NumPy `.npy` and `.npz` files from C++20.

This is a modernised fork of [rogersce/cnpy](https://github.com/rogersce/cnpy). The library was
rewritten for C++20, and the old `npy_load`/`npy_save` API was replaced by a new one (see
[Migrating from cnpy 1.x](#migrating-from-cnpy-1x)).

```cpp
#include <cnpy/cnpy.hpp>

std::vector<double> data(3 * 4);
cnpy::save_npy("a.npy", data, {3, 4});                  // byte-identical to numpy.save

cnpy::ndarray a = cnpy::load_npy("a.npy");
std::span<const double> v = a.data<double>();            // throws cnpy::dtype_error on mismatch
double x = a.at<double>(2, 3);                           // bounds-checked, C or Fortran order

cnpy::ndarray big = cnpy::map_npy("huge.npy");           // memory-mapped, zero copy

cnpy::npz_writer npz("bundle.npz");
npz.add("weights", data, {3, 4}, cnpy::compression::deflate());
npz.add("bias", std::vector<float>{0.1f, 0.2f});
npz.close();

cnpy::npz_reader r("bundle.npz");                        // random access via the zip directory
cnpy::ndarray w = r.load("weights");
```

## Features

- **Paths are `std::filesystem::path`.** They reach the OS as native strings (UTF-16 on Windows), never
  through the C locale, so Unicode file names work everywhere.
- **NumPy-exact output.** `.npy` headers match `numpy.save` byte for byte: 64-byte alignment,
  spare room for growing the array, and format version 1.0 or 2.0. The test suite checks this against
  NumPy.
- **Reads what NumPy writes.** It reads format versions 1.0, 2.0 and 3.0, C and Fortran order,
  big- and little-endian data (converted on load), 0-d and empty arrays, and archives from `savez` and
  `savez_compressed`, including their zip64 records.
- **Data types.** It supports `bool`, all integer widths, `float`, `double`, `long double` and
  `std::complex` of those. `S`, `U` and `V` arrays load as raw bytes. Unsupported types (object,
  structured, datetime) raise a clear `cnpy::format_error`; they are never misread silently.
- **`.npz` support.** Archives are random-access through the zip central directory. Members can be
  deflate-compressed at a chosen level, and CRC-32 is checked on load. zip64 support covers more than
  65,535 members and files larger than 4 GiB. Arrays can be added to an existing archive.
- **Memory mapping** for `.npy` files (`read_only`, `copy_on_write` or `read_write`) and for
  uncompressed `.npz` members.
- **Appending:** `append_npy` grows a file along its growth axis. The header is rewritten in place, and
  the data is written before the header, so an interrupted append leaves the old array intact.
- **Safe types and errors.** Typed `std::span` views check the dtype and alignment. Errors are reported
  as exceptions (`cnpy::io_error`, `format_error`, `dtype_error`). The library never calls `assert`,
  `abort` or prints to `std::cout`.
- **Fast:** no regex, no zero-filled buffers, positional I/O and streaming (de)compression
  (see [Performance](#performance)).

## Building

The library needs a C++20 compiler (GCC 11+, Clang 14+, MSVC 19.30+, Apple Clang 15+), CMake 3.21+ and
zlib. CI tests it on Linux (GCC, Clang), macOS and Windows (MSVC), with ASan and UBSan.

### With Conan 2 (recommended)

```bash
pip install conan numpy            # numpy is only needed for the interop test
conan profile detect
conan install . --build=missing -s compiler.cppstd=20 -s build_type=Release
cmake --preset conan-release       # Visual Studio: cmake --preset conan-default
cmake --build --preset conan-release
ctest --preset conan-release --output-on-failure
```

To use cnpy as a Conan package, run `conan create . -s compiler.cppstd=20`, then add
`self.requires("cnpy/2.0.0")` to your recipe.

### Plain CMake

If zlib is already installed (`zlib1g-dev`, `brew install zlib`, vcpkg, ...):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
cmake --install build --prefix /usr/local
```

Consume it from CMake:

```cmake
find_package(cnpy 2 REQUIRED)
target_link_libraries(my_app PRIVATE cnpy::cnpy)
```

You can also add it with `add_subdirectory(cnpy)` or `FetchContent`.

| CMake option              | Default         | Meaning                                     |
|---------------------------|-----------------|---------------------------------------------|
| `CNPY_BUILD_TESTS`        | top-level only  | doctest unit tests and the NumPy interop test |
| `CNPY_BUILD_EXAMPLES`     | top-level only  | `examples/example.cpp`                      |
| `CNPY_BUILD_BENCHMARKS`   | `OFF`           | `bench/bench.cpp`                           |
| `CNPY_ENABLE_SANITIZERS`  | `OFF`           | AddressSanitizer + UBSan (GCC/Clang)        |
| `CNPY_WARNINGS_AS_ERRORS` | `OFF`           | `-Werror` / `/WX`                           |
| `BUILD_SHARED_LIBS`       | `OFF`           | shared library (symbols exported explicitly) |

## API

Everything is in the `cnpy` namespace and comes from one header, `<cnpy/cnpy.hpp>`.

### `.npy`

| Function | Purpose |
|---|---|
| `ndarray load_npy(path)` | Read the whole file. Foreign byte order is converted to native. |
| `ndarray map_npy(path, map_mode = read_only)` | Memory-map the file. Byte order is kept as stored. |
| `save_npy(path, range, shape, fortran_order = false)` | Write any contiguous range of a supported element type. |
| `save_npy(path, range)` | Write a 1-d array. |
| `save_npy(path, ndarray)` / `save_npy(path, bytes, dtype, shape, fortran)` | Write an array, or raw bytes. |
| `append_npy(path, range, shape, fortran_order = false)` | Append along axis 0 (or along the last axis in Fortran order). Creates the file if it is missing. |

### `.npz`

| API | Purpose |
|---|---|
| `npz_reader r(path, {.memory_map, .verify_crc})` | Open and index an archive. |
| `r.names()`, `r.contains(name)`, `r.size()` | Inspect the archive. |
| `r.load(name)`, `r.load_all()` | Load one array, or all of them as a `std::map`. `load` is thread-safe. |
| `npz_writer w(path, npz_mode::truncate \| append)` | Stream arrays into an archive. |
| `w.add(name, range, shape, compression::deflate(level))` | Add a member. Duplicate names are rejected. |
| `w.close()` | Write the central directory (the destructor calls it too). |
| `load_npz(path)`, `load_npz(path, name)` | Convenience readers. |
| `save_npz(path, name, range, shape, mode, compression)` / `save_npz(path, npz_t)` | Convenience writers. |

When a call has to tell a shape from a compression setting, write the shape explicitly:
`w.add("x", values, cnpy::shape_t{})` for a 0-d array.

### `ndarray`

`dtype()`, `shape()`, `ndim()`, `size()`, `nbytes()`, `fortran_order()`, `is_mapped()`, `writable()`,
`bytes()`, `data<T>()` (read-only `std::span<const T>`), `mutable_data<T>()`, `to_vector<T>()`,
`at<T>(i, j, ...)`, `copy()` (deep copy), `as_native()` (byte-swapped copy if needed),
`ndarray::from(range, shape)`.

Copies of an `ndarray` are shallow, like NumPy views: they share one buffer. Call `copy()` when you need
an independent buffer. Owned buffers are 64-byte aligned.

`std::vector<bool>` is not contiguous, so it cannot be passed directly. Use `std::unique_ptr<bool[]>`,
`std::array<bool, N>` or a `std::span<const bool>` over your own storage instead.

### `dtype`

`dtype::of<T>()`, `dtype::parse("<f8")`, `str()`, `is_native()`, `native()` and `holds<T>()`. Two dtypes
with the same kind and size are compatible: an `int64` array can be read as `long` or `long long`.

## Migrating from cnpy 1.x

| cnpy 1.x | cnpy 2 |
|---|---|
| `#include "cnpy.h"` | `#include <cnpy/cnpy.hpp>` |
| `cnpy::npy_load(std::string)` | `cnpy::load_npy(fs::path)` |
| `cnpy::npy_save(f, ptr, shape, "w")` | `cnpy::save_npy(f, std::span(ptr, n), shape)` |
| `cnpy::npy_save(f, ptr, shape, "a")` | `cnpy::append_npy(f, std::span(ptr, n), shape)` |
| `cnpy::npz_save(z, name, ptr, shape, "w"/"a")` | `cnpy::save_npz(z, name, std::span(ptr, n), shape, npz_mode::truncate/append)` |
| `cnpy::npz_load(z)` / `npz_load(z, name)` | `cnpy::load_npz(z)` / `load_npz(z, name)` |
| `NpyArray::data<T>()` (unchecked `T*`) | `ndarray::data<T>()` (checked span), `mutable_data<T>()` |
| `arr.shape`, `arr.word_size`, `arr.num_vals` | `arr.shape()`, `arr.itemsize()`, `arr.size()` |
| `arr.as_vec<T>()` | `arr.to_vector<T>()` |
| `npz_t = std::map<std::string, NpyArray>` | `npz_t = std::map<std::string, ndarray, std::less<>>` |

The Python 2 helper scripts `mat2npz`, `npy2mat` and `npz2mat` were removed. SciPy's `loadmat` and
`savemat` do the same job in a few lines of Python.

### Bugs fixed along the way

- `npz_load(file, name)` failed with "failed fread" on archives written by `numpy.savez` or
  `numpy.savez_compressed`, because it ignored the zip64 sizes Python stores in the local headers.
  Members are now located through the central directory.
- Big-endian data, out-of-range header fields and corrupt archives were only caught by `assert`. In
  release builds they were silently misread.
- Compressed members were inflated into a temporary buffer and then copied again. They are now streamed
  straight into the array.
- The header shape was parsed with `std::regex` and `stoi`. Shapes beyond `INT_MAX` overflowed.
- CRC-32 values were never verified.
- On Windows, paths went through narrow `fopen`, so non-ASCII file names failed.

## Performance

Run `cnpy_bench [MiB] [dir]` (build with `-DCNPY_BUILD_BENCHMARKS=ON`). The table compares a 256 MiB
`float64` array, page cache warm, on a 2-vCPU Linux VM with GCC 13 and `-O2`. Each figure is the best
of 3 runs, and writes in particular vary a lot from run to run.

| Operation | cnpy 1.x | cnpy 2 |
|---|---:|---:|
| `load_npy` | ~180 ms | ~135 ms |
| `map_npy` + sum over all elements | n/a | ~45 ms |
| `save_npy` | 160–200 ms | 100–170 ms |
| load `.npz` member (stored) | ~185 ms | ~130 ms without CRC check, ~160 ms with |
| load `.npz` member (deflate) | ~365 ms | ~170 ms |
| write `.npz` member (deflate, level 1) | not supported | ~160 ms |

## Testing

- `tests/test_*.cpp` contain the doctest unit tests. They cover every element type, Fortran order,
  byte order, appending (including rewriting legacy headers), memory mapping, CRC and deflate corruption,
  Unicode paths, and zip64 archives with more than 65,535 members.
- `tests/interop/interop.py` runs a NumPy round trip in both directions. NumPy writes about 100 files,
  which C++ reads and checks. C++ then writes files, and NumPy checks the values. Every `.npy` file must
  also be byte-identical to `numpy.save` output, and Python's `zipfile` validates the archives.

## License

MIT. See [LICENSE](LICENSE). Originally written by Carl Rogers (2011).

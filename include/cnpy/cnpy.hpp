// cnpy — read and write NumPy .npy / .npz files from C++20.
//
// Copyright (C) 2011 Carl Rogers, 2026 cnpy contributors
// Released under the MIT License, see LICENSE.

#pragma once

#include <cnpy/export.hpp>

#include <array>
#include <bit>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4251 4275) // std members/bases of exported classes in DLL builds
#endif

namespace cnpy {

namespace fs = std::filesystem;

static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big,
              "cnpy does not support mixed-endian platforms");

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Base class of every exception thrown by cnpy.
class CNPY_EXPORT error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Opening, reading or writing a file failed.
class CNPY_EXPORT io_error : public error {
public:
    io_error(const std::string& what, fs::path path, std::error_code ec = {});
    [[nodiscard]] const fs::path& path() const noexcept { return path_; }
    [[nodiscard]] std::error_code code() const noexcept { return code_; }

private:
    fs::path path_;
    std::error_code code_;
};

/// The file is not a valid .npy/.npz file, or uses a feature cnpy does not support.
class CNPY_EXPORT format_error : public error {
public:
    using error::error;
};

/// An array was accessed as an element type that does not match its dtype.
class CNPY_EXPORT dtype_error : public error {
public:
    using error::error;
};

// ---------------------------------------------------------------------------
// Element types and dtype
// ---------------------------------------------------------------------------

namespace detail {
template <class T>
inline constexpr bool is_complex_v = false;
template <class T>
inline constexpr bool is_complex_v<std::complex<T>> = true;
} // namespace detail

/// Types that map directly onto a NumPy scalar dtype:
/// bool, all integer and character types, float/double/long double and std::complex thereof.
template <class T>
concept element = (std::is_arithmetic_v<T> || detail::is_complex_v<T>) && std::same_as<T, std::remove_cv_t<T>>;

enum class byte_order : char {
    little = '<',
    big = '>',
    not_applicable = '|',
};

inline constexpr byte_order native_byte_order =
    std::endian::native == std::endian::little ? byte_order::little : byte_order::big;

/// A NumPy array-protocol type descriptor such as "<f8", "|b1" or ">c16".
/// Supported kinds: b (bool), i, u, f, c, S (bytes), U (UCS-4 string), V (raw).
struct CNPY_EXPORT dtype {
    char kind = 'V';
    std::size_t itemsize = 0; ///< size of one element in bytes
    byte_order order = byte_order::not_applicable;

    /// dtype of a C++ element type, in native byte order.
    template <element T>
    [[nodiscard]] static constexpr dtype of() noexcept {
        dtype d;
        if constexpr (std::same_as<T, bool>) {
            d.kind = 'b';
        } else if constexpr (std::is_integral_v<T>) {
            d.kind = std::is_signed_v<T> ? 'i' : 'u';
        } else if constexpr (std::is_floating_point_v<T>) {
            d.kind = 'f';
        } else {
            d.kind = 'c';
        }
        d.itemsize = sizeof(T);
        d.order = sizeof(T) == 1 ? byte_order::not_applicable : native_byte_order;
        return d;
    }

    /// Parses a descriptor string ("<f8", "|S10", ">U3", ...). Throws format_error.
    [[nodiscard]] static dtype parse(std::string_view descr);

    /// The descriptor string NumPy would write for this dtype.
    [[nodiscard]] std::string str() const;

    /// True if the data can be read directly on this machine without byte swapping.
    [[nodiscard]] constexpr bool is_native() const noexcept {
        return order == byte_order::not_applicable || order == native_byte_order;
    }

    /// The same type in native byte order.
    [[nodiscard]] constexpr dtype native() const noexcept {
        dtype d = *this;
        if (d.order != byte_order::not_applicable) d.order = native_byte_order;
        return d;
    }

    /// True if an array with this dtype can be viewed as `T` (same kind, size and native order).
    template <element T>
    [[nodiscard]] constexpr bool holds() const noexcept {
        constexpr dtype want = of<T>();
        return kind == want.kind && itemsize == want.itemsize && is_native();
    }

    friend constexpr bool operator==(const dtype&, const dtype&) noexcept = default;
};

using shape_t = std::vector<std::size_t>;

/// Number of elements described by a shape (1 for a 0-d shape). Throws format_error on overflow.
[[nodiscard]] CNPY_EXPORT std::size_t element_count(std::span<const std::size_t> shape);

// ---------------------------------------------------------------------------
// ndarray
// ---------------------------------------------------------------------------

namespace detail {
struct array_access;
}

/// An n-dimensional array loaded from (or to be written to) a .npy file.
///
/// Copies are shallow: they share the same buffer, like NumPy views. Use copy() for a deep copy.
/// The buffer is 64-byte aligned when owned by the array; a memory-mapped array points straight
/// into the mapped file.
class CNPY_EXPORT ndarray {
public:
    ndarray() = default;

    /// Allocates a zero-initialised array.
    ndarray(cnpy::dtype type, shape_t shape, bool fortran_order = false);

    /// Creates an array holding a copy of `values`.
    template <std::ranges::contiguous_range R>
        requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
    [[nodiscard]] static ndarray from(const R& values, shape_t shape, bool fortran_order = false) {
        using T = std::ranges::range_value_t<R>;
        ndarray a(cnpy::dtype::of<T>(), std::move(shape), fortran_order);
        check_size(a.size(), std::ranges::size(values));
        if (a.nbytes() != 0) std::memcpy(a.data_.get(), std::ranges::data(values), a.nbytes());
        return a;
    }

    /// Creates a 1-d array holding a copy of `values`.
    template <std::ranges::contiguous_range R>
        requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
    [[nodiscard]] static ndarray from(const R& values) {
        return from(values, shape_t{std::ranges::size(values)});
    }

    [[nodiscard]] const cnpy::dtype& dtype() const noexcept { return dtype_; }
    [[nodiscard]] const shape_t& shape() const noexcept { return shape_; }
    [[nodiscard]] std::size_t ndim() const noexcept { return shape_.size(); }
    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] std::size_t itemsize() const noexcept { return dtype_.itemsize; }
    [[nodiscard]] std::size_t nbytes() const noexcept { return size_ * dtype_.itemsize; }
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] bool fortran_order() const noexcept { return fortran_order_; }
    [[nodiscard]] bool is_mapped() const noexcept { return mapped_; }
    [[nodiscard]] bool writable() const noexcept { return writable_; }

    [[nodiscard]] std::span<const std::byte> bytes() const noexcept { return {data_.get(), nbytes()}; }
    /// Throws error if the array is a read-only memory mapping.
    [[nodiscard]] std::span<std::byte> mutable_bytes() {
        require_writable();
        return {data_.get(), nbytes()};
    }

    /// Typed read-only view. Throws dtype_error if the dtype does not match `T`.
    template <element T>
    [[nodiscard]] std::span<const T> data() const {
        check_type<T>();
        return {reinterpret_cast<const T*>(data_.get()), size_};
    }

    /// Typed mutable view. Throws dtype_error on type mismatch, error if read-only.
    template <element T>
    [[nodiscard]] std::span<T> mutable_data() {
        check_type<T>();
        require_writable();
        return {reinterpret_cast<T*>(data_.get()), size_};
    }

    template <element T>
    [[nodiscard]] std::vector<T> to_vector() const {
        auto d = data<T>();
        return {d.begin(), d.end()};
    }

    /// Bounds-checked element access honouring C or Fortran order: a.at<float>(i, j, k).
    template <element T, std::integral... I>
    [[nodiscard]] const T& at(I... index) const {
        const std::array<std::size_t, sizeof...(I)> idx{static_cast<std::size_t>(index)...};
        return data<T>()[linear_index(idx)];
    }

    /// Deep copy into a new, owned and writable buffer.
    [[nodiscard]] ndarray copy() const;

    /// Returns *this if the data is already in native byte order, otherwise a byte-swapped copy.
    [[nodiscard]] ndarray as_native() const;

private:
    friend struct detail::array_access;

    template <element T>
    void check_type() const {
        if (!dtype_.holds<T>()) type_mismatch(cnpy::dtype::of<T>());
    }
    [[noreturn]] void type_mismatch(const cnpy::dtype& requested) const;
    void require_writable() const;
    [[nodiscard]] std::size_t linear_index(std::span<const std::size_t> index) const;
    static void check_size(std::size_t expected, std::size_t actual);

    cnpy::dtype dtype_{};
    shape_t shape_{};
    std::size_t size_ = 0;
    bool fortran_order_ = false;
    bool writable_ = true;
    bool mapped_ = false;
    std::shared_ptr<std::byte> data_{};
};

// ---------------------------------------------------------------------------
// .npy
// ---------------------------------------------------------------------------

/// Reads a whole .npy file into memory. Data in foreign byte order is converted to native order.
[[nodiscard]] CNPY_EXPORT ndarray load_npy(const fs::path& path);

enum class map_mode {
    read_only,     ///< mutable access throws
    copy_on_write, ///< writes stay private to this process, the file is untouched
    read_write,    ///< writes go straight to the file
};

/// Memory-maps a .npy file. Pages are loaded lazily by the OS; nothing is copied.
/// The byte order is kept as stored; call as_native() for foreign-endian files.
[[nodiscard]] CNPY_EXPORT ndarray map_npy(const fs::path& path, map_mode mode = map_mode::read_only);

/// Writes raw bytes as a .npy file (the bytes must match dtype and shape).
CNPY_EXPORT void save_npy(const fs::path& path, std::span<const std::byte> bytes, const dtype& type,
                          std::span<const std::size_t> shape, bool fortran_order = false);

CNPY_EXPORT void save_npy(const fs::path& path, const ndarray& array);

/// Writes `values` with the given shape. The output is byte-identical to numpy.save().
template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void save_npy(const fs::path& path, const R& values, const shape_t& shape, bool fortran_order = false) {
    using T = std::ranges::range_value_t<R>;
    const std::span<const T> s(std::ranges::data(values), std::ranges::size(values));
    save_npy(path, std::as_bytes(s), dtype::of<T>(), shape, fortran_order);
}

/// Writes `values` as a 1-d array.
template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void save_npy(const fs::path& path, const R& values) {
    save_npy(path, values, shape_t{std::ranges::size(values)});
}

/// Appends data along the growth axis (axis 0, or the last axis for Fortran-ordered files).
/// Creates the file if it does not exist. The header is rewritten in place whenever possible.
CNPY_EXPORT void append_npy(const fs::path& path, std::span<const std::byte> bytes, const dtype& type,
                            std::span<const std::size_t> shape, bool fortran_order = false);

template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void append_npy(const fs::path& path, const R& values, const shape_t& shape, bool fortran_order = false) {
    using T = std::ranges::range_value_t<R>;
    const std::span<const T> s(std::ranges::data(values), std::ranges::size(values));
    append_npy(path, std::as_bytes(s), dtype::of<T>(), shape, fortran_order);
}

template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void append_npy(const fs::path& path, const R& values) {
    append_npy(path, values, shape_t{std::ranges::size(values)});
}

// ---------------------------------------------------------------------------
// .npz
// ---------------------------------------------------------------------------

using npz_t = std::map<std::string, ndarray, std::less<>>;

struct npz_read_options {
    /// Memory-map the archive. Uncompressed members become zero-copy views when suitably aligned.
    bool memory_map = false;
    /// Verify the CRC-32 of every member that is read into memory.
    bool verify_crc = true;
};

/// Random-access reader for .npz archives (zip, optionally deflate-compressed, zip64 aware).
/// load() is safe to call concurrently from several threads.
class CNPY_EXPORT npz_reader {
public:
    explicit npz_reader(const fs::path& path, npz_read_options options = {});
    ~npz_reader();
    npz_reader(npz_reader&&) noexcept;
    npz_reader& operator=(npz_reader&&) noexcept;
    npz_reader(const npz_reader&) = delete;
    npz_reader& operator=(const npz_reader&) = delete;

    /// Array names in archive order, without the ".npy" suffix.
    [[nodiscard]] std::vector<std::string> names() const;
    [[nodiscard]] bool contains(std::string_view name) const;
    [[nodiscard]] std::size_t size() const noexcept;

    /// Throws error if there is no array called `name`.
    [[nodiscard]] ndarray load(std::string_view name) const;
    [[nodiscard]] npz_t load_all() const;

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

[[nodiscard]] CNPY_EXPORT npz_t load_npz(const fs::path& path, npz_read_options options = {});
[[nodiscard]] CNPY_EXPORT ndarray load_npz(const fs::path& path, std::string_view name, npz_read_options options = {});

/// How a member is stored in the archive.
struct compression {
    int level = 0; ///< 0 = stored (numpy.savez), 1..9 = deflate (numpy.savez_compressed uses 6)

    [[nodiscard]] static constexpr compression none() noexcept { return {0}; }
    [[nodiscard]] static constexpr compression deflate(int level = 6) noexcept { return {level}; }
    friend constexpr bool operator==(compression, compression) noexcept = default;
};

enum class npz_mode {
    truncate, ///< start a new archive
    append,   ///< add to an existing archive (created if missing)
};

/// Streams arrays into a .npz archive. The central directory is written by close() (or the
/// destructor, which swallows errors — call close() explicitly to see them).
class CNPY_EXPORT npz_writer {
public:
    explicit npz_writer(const fs::path& path, npz_mode mode = npz_mode::truncate);
    ~npz_writer();
    npz_writer(npz_writer&&) noexcept;
    npz_writer& operator=(npz_writer&&) noexcept;
    npz_writer(const npz_writer&) = delete;
    npz_writer& operator=(const npz_writer&) = delete;

    /// Adds raw bytes as "<name>.npy". Throws std::invalid_argument for duplicate names.
    void add(std::string_view name, std::span<const std::byte> bytes, const dtype& type,
             std::span<const std::size_t> shape, bool fortran_order = false, compression c = compression::none());

    void add(std::string_view name, const ndarray& array, compression c = compression::none());

    template <std::ranges::contiguous_range R>
        requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
    void add(std::string_view name, const R& values, const shape_t& shape, compression c = compression::none()) {
        using T = std::ranges::range_value_t<R>;
        const std::span<const T> s(std::ranges::data(values), std::ranges::size(values));
        add(name, std::as_bytes(s), dtype::of<T>(), shape, false, c);
    }

    template <std::ranges::contiguous_range R>
        requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
    void add(std::string_view name, const R& values, compression c = compression::none()) {
        add(name, values, shape_t{std::ranges::size(values)}, c);
    }

    /// Writes the central directory and closes the file. Idempotent.
    void close();

private:
    struct impl;
    std::unique_ptr<impl> impl_;
};

/// Writes (or, with npz_mode::append, adds) a single array to an archive.
template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void save_npz(const fs::path& path, std::string_view name, const R& values, const shape_t& shape,
              npz_mode mode = npz_mode::truncate, compression c = compression::none()) {
    npz_writer w(path, mode);
    w.add(name, values, shape, c);
    w.close();
}

template <std::ranges::contiguous_range R>
    requires std::ranges::sized_range<R> && element<std::ranges::range_value_t<R>>
void save_npz(const fs::path& path, std::string_view name, const R& values, npz_mode mode = npz_mode::truncate,
              compression c = compression::none()) {
    save_npz(path, name, values, shape_t{std::ranges::size(values)}, mode, c);
}

/// Writes a whole set of arrays (like numpy.savez / numpy.savez_compressed).
CNPY_EXPORT void save_npz(const fs::path& path, const npz_t& arrays, compression c = compression::none());

} // namespace cnpy

#if defined(_MSC_VER)
#pragma warning(pop)
#endif

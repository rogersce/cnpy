// Internal helpers shared by the cnpy translation units. Not installed.
#pragma once

#include <cnpy/cnpy.hpp>

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <string_view>

namespace cnpy::detail {

/// Path as UTF-8 for error messages (never throws for unrepresentable characters, unlike string()).
[[nodiscard]] inline std::string utf8(const fs::path& p) {
    const auto u = p.u8string();
    return {u.begin(), u.end()};
}

// ----------------------------------------------------------------------------- memory

/// Uninitialised, 64-byte aligned buffer (no zero-fill: we always overwrite it).
[[nodiscard]] std::shared_ptr<std::byte> allocate(std::size_t nbytes);

/// Swaps every element of `bytes` (laid out as `type`) between little and big endian.
void byteswap(std::span<std::byte> bytes, const dtype& type) noexcept;

// ----------------------------------------------------------------------------- files

enum class open_mode { read, read_write, read_write_create, create };

/// Thin RAII wrapper around a native file handle. All I/O is positional (pread/pwrite style),
/// so a const file can be read from several threads at once. Paths go to the OS as native
/// strings (UTF-16 on Windows), never through the C locale.
class file {
public:
    file(const fs::path& path, open_mode mode);
    ~file();
    file(file&& other) noexcept;
    file& operator=(file&& other) noexcept;
    file(const file&) = delete;
    file& operator=(const file&) = delete;

    [[nodiscard]] std::uint64_t size() const;
    /// Reads up to out.size() bytes; returns the number read (short only at end of file).
    std::size_t read_at(std::uint64_t offset, std::span<std::byte> out) const;
    /// Reads exactly out.size() bytes or throws format_error("unexpected end of file").
    void read_exact_at(std::uint64_t offset, std::span<std::byte> out) const;
    void write_at(std::uint64_t offset, std::span<const std::byte> data);
    void truncate(std::uint64_t size);
    void close();

    [[nodiscard]] const fs::path& path() const noexcept { return path_; }
    [[nodiscard]] std::intptr_t native_handle() const noexcept { return handle_; }

private:
    fs::path path_;
    std::intptr_t handle_ = -1;
};

/// Maps a whole file and returns a pointer to its first byte that keeps the mapping alive.
[[nodiscard]] std::shared_ptr<std::byte> map_file(const file& f, map_mode mode, std::uint64_t size);

[[noreturn]] void throw_io(const std::string& what, const fs::path& path, std::error_code ec);
[[noreturn]] void throw_last_io_error(const std::string& what, const fs::path& path);

// ----------------------------------------------------------------------------- npy header

inline constexpr std::size_t npy_prefix_v1 = 10; // magic(6) + version(2) + uint16 length
inline constexpr std::size_t npy_prefix_v2 = 12; // magic(6) + version(2) + uint32 length
inline constexpr std::size_t npy_align = 64;
inline constexpr std::size_t growth_axis_max_digits = 21;

struct header_info {
    dtype type;
    shape_t shape;
    bool fortran_order = false;
    std::uint8_t major_version = 1;
    std::size_t header_size = 0; ///< prefix + dictionary, i.e. offset of the data
    std::size_t nbytes = 0;      ///< size of the array data
};

/// The Python dict literal (without padding), e.g. "{'descr': '<f8', 'fortran_order': False, 'shape': (3,), }".
[[nodiscard]] std::string header_dict(const dtype& type, std::span<const std::size_t> shape, bool fortran_order);

/// A complete header exactly as numpy.save would write it (growth padding, 64-byte alignment, v1/v2).
[[nodiscard]] std::string make_header(const dtype& type, std::span<const std::size_t> shape, bool fortran_order);

/// A header of exactly `total_size` bytes and the given major version, or an empty string if the
/// dictionary does not fit. Used to rewrite headers of existing files in place.
[[nodiscard]] std::string make_header_of_size(const dtype& type, std::span<const std::size_t> shape, bool fortran_order,
                                              std::uint8_t major_version, std::size_t total_size);

/// Length of the fixed prefix and of the full header given the first 12 bytes (or fewer if the
/// file is smaller). Throws format_error on bad magic/version.
struct header_prefix {
    std::uint8_t major_version;
    std::size_t prefix_size;
    std::size_t header_size;
};
[[nodiscard]] header_prefix parse_header_prefix(std::span<const std::byte> first_bytes);

/// Parses the dictionary part of a header (between prefix and data).
[[nodiscard]] header_info parse_header_dict(std::string_view dict);

/// Reads and parses the header of a .npy file stored at `offset` in `f`, limited to `available`
/// bytes (the file size, or the member size inside a zip).
[[nodiscard]] header_info read_header(const file& f, std::uint64_t offset, std::uint64_t available);

/// Parses a header from memory.
[[nodiscard]] header_info read_header(std::span<const std::byte> bytes);

/// Throws std::invalid_argument unless bytes.size() == element_count(shape) * type.itemsize.
void check_bytes(std::span<const std::byte> bytes, const dtype& type, std::span<const std::size_t> shape);

// ----------------------------------------------------------------------------- ndarray internals

struct array_access {
    static ndarray make(const dtype& type, shape_t shape, bool fortran_order, std::shared_ptr<std::byte> data,
                        bool writable, bool mapped);
    static std::byte* raw(ndarray& a) noexcept { return a.data_.get(); }
};

/// Allocates an uninitialised owned array (to be filled by the caller).
[[nodiscard]] ndarray make_uninitialized(const header_info& h);

/// Little-endian integer helpers for the zip structures.
template <std::unsigned_integral T>
[[nodiscard]] inline T load_le(const std::byte* p) noexcept {
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < sizeof(T); ++i)
        v |= std::uint64_t{std::to_integer<std::uint8_t>(p[i])} << (8 * i);
    return static_cast<T>(v);
}

template <std::unsigned_integral T>
inline void store_le(std::string& out, T v) {
    for (std::size_t i = 0; i < sizeof(T); ++i)
        out.push_back(static_cast<char>((v >> (8 * i)) & 0xFF));
}

[[nodiscard]] inline std::span<const std::byte> as_bytes(std::string_view s) noexcept {
    return std::as_bytes(std::span<const char>(s.data(), s.size()));
}

} // namespace cnpy::detail

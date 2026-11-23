// .npy header handling, load/map/save/append.
#include "detail.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <limits>

namespace cnpy {

namespace detail {

namespace {

constexpr std::string_view magic = "\x93NUMPY";

// --------------------------------------------------------------------------- writing headers

std::string shape_repr(std::span<const std::size_t> shape) {
    // Mirrors Python's repr(tuple): (), (3,), (3, 4)
    std::string s = "(";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i != 0) s += ", ";
        s += std::to_string(shape[i]);
    }
    if (shape.size() == 1) s += ',';
    s += ')';
    return s;
}

std::string wrap(std::string_view dict, std::uint8_t major, std::size_t dict_and_padding) {
    // dict_and_padding counts dict + spaces + trailing newline
    std::string out;
    out.reserve(npy_prefix_v2 + dict_and_padding);
    out += magic;
    out.push_back(static_cast<char>(major));
    out.push_back('\0');
    if (major == 1) {
        store_le(out, static_cast<std::uint16_t>(dict_and_padding));
    } else {
        store_le(out, static_cast<std::uint32_t>(dict_and_padding));
    }
    out += dict;
    out.append(dict_and_padding - dict.size() - 1, ' ');
    out.push_back('\n');
    return out;
}

// --------------------------------------------------------------------------- parsing headers

/// Minimal parser for the Python dict literal in a .npy header. No regex, no allocation per token.
class literal_parser {
public:
    explicit literal_parser(std::string_view s) : s_(s) {}

    void skip_ws() noexcept {
        while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\n' || s_[i_] == '\r'))
            ++i_;
    }
    [[nodiscard]] bool at_end() noexcept {
        skip_ws();
        return i_ == s_.size();
    }
    [[nodiscard]] char peek() noexcept {
        skip_ws();
        return i_ < s_.size() ? s_[i_] : '\0';
    }
    bool consume(char c) noexcept {
        if (peek() != c) return false;
        ++i_;
        return true;
    }
    void expect(char c) {
        if (!consume(c)) fail(std::string("expected '") + c + "'");
    }

    std::string_view string() {
        const char quote = peek();
        if (quote != '\'' && quote != '"') fail("expected a string");
        const std::size_t start = ++i_;
        while (i_ < s_.size() && s_[i_] != quote) {
            if (s_[i_] == '\\') ++i_;
            ++i_;
        }
        if (i_ >= s_.size()) fail("unterminated string");
        return s_.substr(start, i_++ - start);
    }

    bool boolean() {
        skip_ws();
        if (s_.substr(i_, 4) == "True") {
            i_ += 4;
            return true;
        }
        if (s_.substr(i_, 5) == "False") {
            i_ += 5;
            return false;
        }
        fail("expected True or False");
    }

    shape_t tuple() {
        expect('(');
        shape_t shape;
        while (!consume(')')) {
            skip_ws();
            std::size_t v = 0;
            const auto [ptr, ec] = std::from_chars(s_.data() + i_, s_.data() + s_.size(), v);
            if (ec != std::errc{}) fail("expected a non-negative integer in shape");
            i_ = static_cast<std::size_t>(ptr - s_.data());
            if (i_ < s_.size() && s_[i_] == 'L') ++i_; // Python 2 long suffix
            shape.push_back(v);
            if (!consume(',')) {
                expect(')');
                break;
            }
        }
        return shape;
    }

    [[noreturn]] void fail(const std::string& why) const {
        throw format_error("malformed .npy header (" + why + " at offset " + std::to_string(i_) +
                           "): " + std::string(s_.substr(0, 200)));
    }

private:
    std::string_view s_;
    std::size_t i_ = 0;
};

std::size_t header_nbytes(const dtype& type, std::span<const std::size_t> shape) {
    const std::size_t n = element_count(shape);
    if (type.itemsize != 0 && n > std::numeric_limits<std::size_t>::max() / type.itemsize)
        throw format_error("array byte size overflows size_t");
    return n * type.itemsize;
}

} // namespace

std::string header_dict(const dtype& type, std::span<const std::size_t> shape, bool fortran_order) {
    // Keys in sorted order, formatted exactly like numpy.lib.format._write_array_header.
    std::string d = "{'descr': '";
    d += type.str();
    d += "', 'fortran_order': ";
    d += fortran_order ? "True" : "False";
    d += ", 'shape': ";
    d += shape_repr(shape);
    d += ", }";
    return d;
}

std::string make_header(const dtype& type, std::span<const std::size_t> shape, bool fortran_order) {
    std::string dict = header_dict(type, shape, fortran_order);
    // Spare room so the growth axis can be enlarged in place (numpy >= 1.24 does the same).
    if (!shape.empty()) {
        const std::size_t digits = std::to_string(fortran_order ? shape.back() : shape.front()).size();
        if (digits < growth_axis_max_digits) dict.append(growth_axis_max_digits - digits, ' ');
    }
    const std::size_t hlen = dict.size() + 1; // + '\n'
    for (const std::uint8_t major : {std::uint8_t{1}, std::uint8_t{2}}) {
        const std::size_t prefix = major == 1 ? npy_prefix_v1 : npy_prefix_v2;
        const std::size_t padlen = npy_align - ((prefix + hlen) % npy_align); // 1..64, as numpy does
        const std::size_t body = hlen + padlen;
        if (major == 1 && body > std::numeric_limits<std::uint16_t>::max()) continue;
        if (body > std::numeric_limits<std::uint32_t>::max()) break;
        return wrap(dict, major, body);
    }
    throw format_error("array header too large");
}

std::string make_header_of_size(const dtype& type, std::span<const std::size_t> shape, bool fortran_order,
                                std::uint8_t major_version, std::size_t total_size) {
    const std::string dict = header_dict(type, shape, fortran_order);
    const std::size_t prefix = major_version == 1 ? npy_prefix_v1 : npy_prefix_v2;
    if (total_size < prefix + dict.size() + 1) return {};
    return wrap(dict, major_version, total_size - prefix);
}

header_prefix parse_header_prefix(std::span<const std::byte> b) {
    if (b.size() < npy_prefix_v1 || std::memcmp(b.data(), magic.data(), magic.size()) != 0)
        throw format_error("not a .npy file (bad magic string)");
    const auto major = std::to_integer<std::uint8_t>(b[6]);
    if (major == 1) {
        const auto len = load_le<std::uint16_t>(b.data() + 8);
        return {major, npy_prefix_v1, npy_prefix_v1 + len};
    }
    if (major == 2 || major == 3) {
        if (b.size() < npy_prefix_v2) throw format_error("truncated .npy header");
        const auto len = load_le<std::uint32_t>(b.data() + 8);
        return {major, npy_prefix_v2, npy_prefix_v2 + std::size_t{len}};
    }
    throw format_error("unsupported .npy format version " + std::to_string(major));
}

header_info parse_header_dict(std::string_view dict) {
    literal_parser p(dict);
    header_info h;
    bool have_descr = false;
    bool have_order = false;
    bool have_shape = false;

    p.expect('{');
    while (!p.consume('}')) {
        const std::string_view key = p.string();
        p.expect(':');
        if (key == "descr") {
            if (p.peek() == '[') throw format_error("structured dtypes (record arrays) are not supported");
            h.type = dtype::parse(p.string());
            have_descr = true;
        } else if (key == "fortran_order") {
            h.fortran_order = p.boolean();
            have_order = true;
        } else if (key == "shape") {
            h.shape = p.tuple();
            have_shape = true;
        } else {
            p.fail("unexpected key '" + std::string(key) + "'");
        }
        if (!p.consume(',')) {
            p.expect('}');
            break;
        }
    }
    if (!p.at_end()) p.fail("trailing characters");
    if (!have_descr || !have_order || !have_shape) p.fail("missing 'descr', 'fortran_order' or 'shape'");
    h.nbytes = header_nbytes(h.type, h.shape);
    return h;
}

header_info read_header(std::span<const std::byte> bytes) {
    const header_prefix pre = parse_header_prefix(bytes.first(std::min(bytes.size(), npy_prefix_v2)));
    if (pre.header_size > bytes.size()) throw format_error("truncated .npy header");
    const auto dict = bytes.subspan(pre.prefix_size, pre.header_size - pre.prefix_size);
    header_info h = parse_header_dict(std::string_view(reinterpret_cast<const char*>(dict.data()), dict.size()));
    h.major_version = pre.major_version;
    h.header_size = pre.header_size;
    return h;
}

header_info read_header(const file& f, std::uint64_t offset, std::uint64_t available) {
    std::array<std::byte, npy_prefix_v2> first{};
    const std::size_t want = static_cast<std::size_t>(std::min<std::uint64_t>(available, first.size()));
    const std::size_t got = f.read_at(offset, std::span(first).first(want));
    const header_prefix pre = parse_header_prefix(std::span(first).first(got));
    if (pre.header_size > available) throw format_error("truncated .npy header in '" + utf8(f.path()) + "'");

    std::string dict(pre.header_size - pre.prefix_size, '\0');
    f.read_exact_at(offset + pre.prefix_size, std::as_writable_bytes(std::span(dict)));
    header_info h = parse_header_dict(dict);
    h.major_version = pre.major_version;
    h.header_size = pre.header_size;
    return h;
}

void check_bytes(std::span<const std::byte> bytes, const dtype& type, std::span<const std::size_t> shape) {
    const std::size_t n = element_count(shape);
    if (type.itemsize != 0 && n > std::numeric_limits<std::size_t>::max() / type.itemsize)
        throw std::invalid_argument("array byte size overflows size_t");
    if (bytes.size() != n * type.itemsize)
        throw std::invalid_argument("data size (" + std::to_string(bytes.size()) +
                                    " bytes) does not match shape and dtype (" + std::to_string(n) + " x " +
                                    std::to_string(type.itemsize) + " bytes)");
}

} // namespace detail

namespace {

using detail::file;
using detail::header_info;
using detail::open_mode;

void check_data_size(const file& f, const header_info& h, std::uint64_t file_size) {
    if (file_size - h.header_size < h.nbytes)
        throw format_error("'" + detail::utf8(f.path()) + "' is truncated: expected " + std::to_string(h.nbytes) +
                           " data bytes, found " + std::to_string(file_size - h.header_size));
}

void write_npy(const fs::path& path, std::span<const std::byte> bytes, const dtype& type,
               std::span<const std::size_t> shape, bool fortran_order) {
    const std::string header = detail::make_header(type, shape, fortran_order);
    file f(path, open_mode::create);
    f.write_at(0, detail::as_bytes(header));
    f.write_at(header.size(), bytes);
}

} // namespace

ndarray load_npy(const fs::path& path) {
    const file f(path, open_mode::read);
    const std::uint64_t size = f.size();
    const header_info h = detail::read_header(f, 0, size);
    check_data_size(f, h, size);

    auto buffer = detail::allocate(h.nbytes);
    f.read_exact_at(h.header_size, {buffer.get(), h.nbytes});
    dtype type = h.type;
    if (!type.is_native()) {
        detail::byteswap({buffer.get(), h.nbytes}, type);
        type = type.native();
    }
    return detail::array_access::make(type, h.shape, h.fortran_order, std::move(buffer), true, false);
}

ndarray map_npy(const fs::path& path, map_mode mode) {
    const file f(path, mode == map_mode::read_write ? open_mode::read_write : open_mode::read);
    const std::uint64_t size = f.size();
    const header_info h = detail::read_header(f, 0, size);
    check_data_size(f, h, size);

    const std::shared_ptr<std::byte> base = detail::map_file(f, mode, size);
    std::shared_ptr<std::byte> data(base, base.get() + h.header_size); // aliasing: keeps mapping alive
    return detail::array_access::make(h.type, h.shape, h.fortran_order, std::move(data), mode != map_mode::read_only,
                                      true);
}

void save_npy(const fs::path& path, std::span<const std::byte> bytes, const dtype& type,
              std::span<const std::size_t> shape, bool fortran_order) {
    detail::check_bytes(bytes, type, shape);
    write_npy(path, bytes, type, shape, fortran_order);
}

void save_npy(const fs::path& path, const ndarray& array) {
    write_npy(path, array.bytes(), array.dtype(), array.shape(), array.fortran_order());
}

void append_npy(const fs::path& path, std::span<const std::byte> bytes, const dtype& type,
                std::span<const std::size_t> shape, bool fortran_order) {
    detail::check_bytes(bytes, type, shape);
    if (shape.empty()) throw std::invalid_argument("append_npy: cannot append a 0-d array");

    std::error_code ec;
    if (!fs::exists(path, ec)) {
        if (ec) detail::throw_io("cannot access file", path, ec);
        write_npy(path, bytes, type, shape, fortran_order);
        return;
    }

    file f(path, open_mode::read_write);
    const std::uint64_t size = f.size();
    const header_info h = detail::read_header(f, 0, size);

    if (h.type != type)
        throw dtype_error("append_npy: file has dtype '" + h.type.str() + "', appended data is '" + type.str() + "'");
    if (h.fortran_order != fortran_order) throw error("append_npy: memory order (C/Fortran) differs from the file");
    if (h.shape.size() != shape.size())
        throw error("append_npy: file has " + std::to_string(h.shape.size()) + " dimensions, appended data has " +
                    std::to_string(shape.size()));
    const std::size_t growth_axis = fortran_order ? shape.size() - 1 : 0;
    for (std::size_t d = 0; d < shape.size(); ++d)
        if (d != growth_axis && h.shape[d] != shape[d])
            throw error("append_npy: shape mismatch on axis " + std::to_string(d));
    if (size != h.header_size + h.nbytes) throw format_error("append_npy: unexpected trailing data in file");

    shape_t new_shape = h.shape;
    new_shape[growth_axis] += shape[growth_axis];
    if (new_shape[growth_axis] < shape[growth_axis]) throw format_error("append_npy: shape overflow");
    (void)element_count(new_shape);

    std::string header = detail::make_header(type, new_shape, fortran_order);
    if (header.size() != h.header_size)
        header = detail::make_header_of_size(type, new_shape, fortran_order, h.major_version, h.header_size);

    if (!header.empty()) {
        // Data first, header last: an interrupted append leaves the old array intact.
        f.write_at(size, bytes);
        f.write_at(0, detail::as_bytes(header));
        return;
    }

    // Legacy file without spare header room: rewrite through a temporary file.
    header = detail::make_header(type, new_shape, fortran_order);
    fs::path tmp = path;
    tmp += ".cnpy-tmp";
    try {
        file out(tmp, open_mode::create);
        out.write_at(0, detail::as_bytes(header));
        std::uint64_t pos = header.size();
        const std::size_t chunk_size = std::size_t{16} << 20;
        auto chunk = detail::allocate(chunk_size);
        for (std::uint64_t off = h.header_size; off < size;) {
            const auto n = static_cast<std::size_t>(std::min<std::uint64_t>(chunk_size, size - off));
            f.read_exact_at(off, {chunk.get(), n});
            out.write_at(pos, {chunk.get(), n});
            off += n;
            pos += n;
        }
        out.write_at(pos, bytes);
        out.close();
        f.close();
        fs::rename(tmp, path, ec);
        if (ec) detail::throw_io("cannot replace file", path, ec);
    } catch (...) {
        fs::remove(tmp, ec);
        throw;
    }
}

} // namespace cnpy

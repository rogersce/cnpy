// dtype and ndarray.
#include "detail.hpp"

#include <algorithm>
#include <charconv>
#include <limits>

namespace cnpy {

// ----------------------------------------------------------------------------- dtype

namespace {

bool valid_itemsize(char kind, std::size_t n) noexcept {
    switch (kind) {
    case 'b': return n == 1;
    case 'i':
    case 'u': return n == 1 || n == 2 || n == 4 || n == 8;
    case 'f': return n == 2 || n == 4 || n == 8 || n == 10 || n == 12 || n == 16;
    case 'c': return n == 8 || n == 16 || n == 20 || n == 24 || n == 32;
    case 'S':
    case 'V': return true;
    case 'U': return n % 4 == 0;
    default: return false;
    }
}

} // namespace

dtype dtype::parse(std::string_view descr) {
    const std::string original(descr);
    auto fail = [&](const char* why) -> format_error {
        return format_error("unsupported dtype descriptor '" + original + "': " + why);
    };

    dtype d;
    byte_order order = native_byte_order;
    bool explicit_order = false;
    if (!descr.empty() &&
        (descr.front() == '<' || descr.front() == '>' || descr.front() == '|' || descr.front() == '=')) {
        explicit_order = true;
        switch (descr.front()) {
        case '<': order = byte_order::little; break;
        case '>': order = byte_order::big; break;
        case '|': order = byte_order::not_applicable; break;
        default: order = native_byte_order; break;
        }
        descr.remove_prefix(1);
    }
    if (descr.empty()) throw fail("missing type code");

    char kind = descr.front();
    if (kind == '?') kind = 'b';
    if (kind == 'a') kind = 'S';
    if (kind == 'O') throw fail("object arrays (pickled Python objects) cannot be read");
    if (kind == 'M' || kind == 'm') throw fail("datetime64/timedelta64 are not supported");
    descr.remove_prefix(1);

    std::size_t count = 0;
    if (descr.empty()) {
        if (kind != 'b') throw fail("missing item size");
        count = 1;
    } else {
        const auto* first = descr.data();
        const auto* last = descr.data() + descr.size();
        const auto [ptr, ec] = std::from_chars(first, last, count);
        if (ec != std::errc{} || ptr != last) throw fail("malformed item size");
    }

    std::size_t itemsize = count;
    if (kind == 'U') {
        if (count > std::numeric_limits<std::size_t>::max() / 4) throw fail("item size too large");
        itemsize = count * 4;
    }
    if (!valid_itemsize(kind, itemsize)) throw fail("unknown type code or item size");

    d.kind = kind;
    d.itemsize = itemsize;
    const bool order_irrelevant =
        kind == 'S' || kind == 'V' || ((kind == 'b' || kind == 'i' || kind == 'u') && itemsize == 1);
    if (order_irrelevant) {
        d.order = byte_order::not_applicable;
    } else if (!explicit_order || order == byte_order::not_applicable) {
        d.order = native_byte_order; // numpy treats "|f8" or "f8" as native
    } else {
        d.order = order;
    }
    return d;
}

std::string dtype::str() const {
    std::string s;
    s.push_back(static_cast<char>(order));
    s.push_back(kind);
    s += std::to_string(kind == 'U' ? itemsize / 4 : itemsize);
    return s;
}

std::size_t element_count(std::span<const std::size_t> shape) {
    std::size_t n = 1;
    for (const std::size_t dim : shape) {
        if (dim != 0 && n > std::numeric_limits<std::size_t>::max() / dim)
            throw format_error("array shape overflows size_t");
        n *= dim;
    }
    return n;
}

namespace detail {

namespace {
std::size_t checked_nbytes(const dtype& type, std::span<const std::size_t> shape) {
    const std::size_t n = element_count(shape);
    if (type.itemsize != 0 && n > std::numeric_limits<std::size_t>::max() / type.itemsize)
        throw format_error("array byte size overflows size_t");
    return n * type.itemsize;
}

template <std::size_t N>
void swap_fixed(std::byte* p, std::size_t count) noexcept {
    for (std::size_t i = 0; i < count; ++i, p += N)
        std::reverse(p, p + N);
}
} // namespace

void byteswap(std::span<std::byte> bytes, const dtype& type) noexcept {
    std::size_t unit = type.itemsize;
    if (type.kind == 'c') unit = type.itemsize / 2; // real and imaginary part separately
    if (type.kind == 'U') unit = 4;
    if (type.kind == 'S' || type.kind == 'V' || type.kind == 'b' || unit <= 1) return;

    std::byte* p = bytes.data();
    const std::size_t count = bytes.size() / unit;
    // Fixed-size reverses compile down to bswap instructions.
    switch (unit) {
    case 2: swap_fixed<2>(p, count); break;
    case 4: swap_fixed<4>(p, count); break;
    case 8: swap_fixed<8>(p, count); break;
    case 16: swap_fixed<16>(p, count); break;
    default:
        for (std::size_t i = 0; i < count; ++i, p += unit)
            std::reverse(p, p + unit);
        break;
    }
}

ndarray array_access::make(const dtype& type, shape_t shape, bool fortran_order, std::shared_ptr<std::byte> data,
                           bool writable, bool mapped) {
    ndarray a;
    a.size_ = element_count(shape);
    (void)checked_nbytes(type, shape);
    a.dtype_ = type;
    a.shape_ = std::move(shape);
    a.fortran_order_ = fortran_order;
    a.data_ = std::move(data);
    a.writable_ = writable;
    a.mapped_ = mapped;
    return a;
}

ndarray make_uninitialized(const header_info& h) {
    return array_access::make(h.type, h.shape, h.fortran_order, allocate(h.nbytes), true, false);
}

} // namespace detail

// ----------------------------------------------------------------------------- ndarray

ndarray::ndarray(cnpy::dtype type, shape_t shape, bool fortran_order)
    : dtype_(type), shape_(std::move(shape)), size_(element_count(shape_)), fortran_order_(fortran_order) {
    const std::size_t n = detail::checked_nbytes(dtype_, shape_);
    data_ = detail::allocate(n);
    if (n != 0) std::memset(data_.get(), 0, n);
}

ndarray ndarray::copy() const {
    ndarray a = detail::array_access::make(dtype_, shape_, fortran_order_, detail::allocate(nbytes()), true, false);
    if (nbytes() != 0) std::memcpy(a.data_.get(), data_.get(), nbytes());
    return a;
}

ndarray ndarray::as_native() const {
    if (dtype_.is_native()) return *this;
    ndarray a = copy();
    detail::byteswap(a.mutable_bytes(), a.dtype_);
    a.dtype_ = a.dtype_.native();
    return a;
}

void ndarray::type_mismatch(const cnpy::dtype& requested) const {
    std::string msg = "array has dtype '" + dtype_.str() + "' but was accessed as '" + requested.str() + "'";
    if (!dtype_.is_native()) msg += " (foreign byte order: use as_native())";
    throw dtype_error(msg);
}

void ndarray::require_writable() const {
    if (!writable_) throw error("array is a read-only memory mapping");
}

void ndarray::check_size(std::size_t expected, std::size_t actual) {
    if (expected != actual)
        throw std::invalid_argument("shape describes " + std::to_string(expected) + " elements but " +
                                    std::to_string(actual) + " were given");
}

std::size_t ndarray::linear_index(std::span<const std::size_t> index) const {
    if (index.size() != shape_.size())
        throw std::out_of_range("array has " + std::to_string(shape_.size()) + " dimensions, got " +
                                std::to_string(index.size()) + " indices");
    for (std::size_t d = 0; d < index.size(); ++d)
        if (index[d] >= shape_[d])
            throw std::out_of_range("index " + std::to_string(index[d]) + " out of range for axis " +
                                    std::to_string(d) + " with size " + std::to_string(shape_[d]));
    std::size_t linear = 0;
    if (fortran_order_) {
        for (std::size_t d = index.size(); d-- > 0;)
            linear = linear * shape_[d] + index[d];
    } else {
        for (std::size_t d = 0; d < index.size(); ++d)
            linear = linear * shape_[d] + index[d];
    }
    return linear;
}

} // namespace cnpy

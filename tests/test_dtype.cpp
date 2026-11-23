#include <cnpy/cnpy.hpp>

#include <doctest/doctest.h>

#include <complex>
#include <cstdint>

using cnpy::byte_order;
using cnpy::dtype;

namespace {
constexpr char native = static_cast<char>(cnpy::native_byte_order);
std::string n(const char* rest) {
    return std::string(1, native) + rest;
}
} // namespace

TEST_CASE("dtype::of maps C++ types to NumPy descriptors") {
    CHECK(dtype::of<bool>().str() == "|b1");
    CHECK(dtype::of<std::int8_t>().str() == "|i1");
    CHECK(dtype::of<std::uint8_t>().str() == "|u1");
    CHECK(dtype::of<std::int16_t>().str() == n("i2"));
    CHECK(dtype::of<std::uint32_t>().str() == n("u4"));
    CHECK(dtype::of<std::int64_t>().str() == n("i8"));
    CHECK(dtype::of<float>().str() == n("f4"));
    CHECK(dtype::of<double>().str() == n("f8"));
    CHECK(dtype::of<std::complex<float>>().str() == n("c8"));
    CHECK(dtype::of<std::complex<double>>().str() == n("c16"));
    CHECK(dtype::of<char>().kind == (std::is_signed_v<char> ? 'i' : 'u'));

    static_assert(dtype::of<double>().holds<double>());
    static_assert(!dtype::of<double>().holds<float>());
    static_assert(dtype::of<long long>().holds<std::int64_t>()); // same kind and size
}

TEST_CASE("dtype::parse round-trips NumPy descriptors") {
    for (const char* s : {"<f8", ">f4", "|b1", "|i1", "|u1", "<i2", ">u8", "<c16", ">c8", "|S10", "<U3", ">U1", "|V16",
                          "<f2", "<f16"}) {
        CAPTURE(s);
        CHECK(dtype::parse(s).str() == s);
    }
}

TEST_CASE("dtype::parse normalises byte order") {
    CHECK(dtype::parse("=f8") == dtype::of<double>());
    CHECK(dtype::parse("|f8") == dtype::of<double>());
    CHECK(dtype::parse("f8") == dtype::of<double>());
    CHECK(dtype::parse("<i1").order == byte_order::not_applicable);
    CHECK(dtype::parse(">S4").order == byte_order::not_applicable);
    CHECK(dtype::parse("?").kind == 'b');
    CHECK(dtype::parse("<U5").itemsize == 20);
}

TEST_CASE("dtype::parse rejects unsupported descriptors") {
    for (const char* s : {"", "<", "<f", "<i3", "<x4", "|O8", "<M8[ns]", "<m8", "<f8x", "<U3a", "<c4"}) {
        CAPTURE(s);
        CHECK_THROWS_AS((void)dtype::parse(s), cnpy::format_error);
    }
}

TEST_CASE("dtype native() and is_native()") {
    const dtype be = dtype::parse(">i4");
    const dtype le = dtype::parse("<i4");
    CHECK(be.is_native() == (cnpy::native_byte_order == byte_order::big));
    CHECK(le.is_native() == (cnpy::native_byte_order == byte_order::little));
    CHECK(be.native() == dtype::of<std::int32_t>());
    CHECK(dtype::parse("|S3").is_native());
}

TEST_CASE("element_count") {
    CHECK(cnpy::element_count(cnpy::shape_t{}) == 1);
    CHECK(cnpy::element_count(cnpy::shape_t{3, 0, 2}) == 0);
    CHECK(cnpy::element_count(cnpy::shape_t{3, 4, 5}) == 60);
    CHECK_THROWS_AS((void)cnpy::element_count(cnpy::shape_t{std::size_t(1) << 40, std::size_t(1) << 40}),
                    cnpy::format_error);
}

#include "helpers.hpp"

#include <doctest/doctest.h>

#include <complex>
#include <cstdint>

using test::temp_dir;
using test::values;

#define CNPY_ELEMENT_TYPES                                                                                             \
    bool, char, signed char, unsigned char, std::int16_t, std::uint16_t, std::int32_t, std::uint32_t, std::int64_t,    \
        std::uint64_t, long, unsigned long long, float, double, long double, std::complex<float>, std::complex<double>

TEST_CASE_TEMPLATE("save_npy / load_npy round-trip", T, CNPY_ELEMENT_TYPES) {
    temp_dir dir;
    const auto path = dir / "a.npy";
    const auto v = values<T>(2 * 3 * 4);
    cnpy::save_npy(path, v, {2, 3, 4});

    const cnpy::ndarray a = cnpy::load_npy(path);
    CHECK(a.dtype() == cnpy::dtype::of<T>());
    CHECK(a.shape() == cnpy::shape_t{2, 3, 4});
    CHECK(a.ndim() == 3);
    CHECK(a.size() == 24);
    CHECK(a.nbytes() == 24 * sizeof(T));
    CHECK_FALSE(a.fortran_order());
    CHECK_FALSE(a.is_mapped());
    CHECK(test::same(a.data<T>(), v));
    CHECK(a.at<T>(1, 2, 3) == v[23]);
    CHECK(a.at<T>(1, 0, 2) == v[1 * 12 + 0 * 4 + 2]);
}

TEST_CASE("header layout follows the NumPy format specification") {
    temp_dir dir;
    const auto path = dir / "h.npy";
    cnpy::save_npy(path, std::vector<double>(12), {3, 4});
    const std::string f = test::read_file(path);

    REQUIRE(f.size() > 10);
    CHECK(f.substr(0, 6) == "\x93NUMPY");
    CHECK(f[6] == 1);
    CHECK(f[7] == 0);
    const std::size_t hlen = test::le16(f, 8);
    CHECK((10 + hlen) % 64 == 0);
    CHECK(f.size() == 10 + hlen + 12 * sizeof(double));
    const std::string dict = f.substr(10, hlen);
    CHECK(dict.rfind("{'descr': '<f8', 'fortran_order': False, 'shape': (3, 4), }", 0) == 0);
    CHECK(dict.back() == '\n');
}

TEST_CASE("0-d, 1-d and empty arrays") {
    temp_dir dir;
    SUBCASE("scalar") {
        const double x = 4.5;
        cnpy::save_npy(dir / "s.npy", std::span(&x, 1), {});
        const auto a = cnpy::load_npy(dir / "s.npy");
        CHECK(a.shape().empty());
        CHECK(a.size() == 1);
        CHECK(a.at<double>() == 4.5);
        CHECK(test::read_file(dir / "s.npy").find("'shape': ()") != std::string::npos);
    }
    SUBCASE("1-d convenience overload") {
        const std::vector<int> v{1, 2, 3};
        cnpy::save_npy(dir / "v.npy", v);
        const auto a = cnpy::load_npy(dir / "v.npy");
        CHECK(a.shape() == cnpy::shape_t{3});
        CHECK(test::read_file(dir / "v.npy").find("'shape': (3,)") != std::string::npos);
    }
    SUBCASE("empty") {
        cnpy::save_npy(dir / "e.npy", std::vector<float>{}, {0, 5});
        const auto a = cnpy::load_npy(dir / "e.npy");
        CHECK(a.shape() == cnpy::shape_t{0, 5});
        CHECK(a.empty());
        CHECK(a.data<float>().empty());
    }
    SUBCASE("std::array and C arrays work as ranges") {
        const std::array<std::uint16_t, 4> arr{1, 2, 3, 4};
        cnpy::save_npy(dir / "arr.npy", arr, {2, 2});
        const float c[3] = {1.f, 2.f, 3.f};
        cnpy::save_npy(dir / "c.npy", c);
        CHECK(cnpy::load_npy(dir / "arr.npy").at<std::uint16_t>(1, 0) == 3);
        CHECK(cnpy::load_npy(dir / "c.npy").to_vector<float>() == std::vector<float>{1.f, 2.f, 3.f});
    }
}

TEST_CASE("Fortran order is stored and honoured by at()") {
    temp_dir dir;
    // Column-major 2x3 matrix [[0, 1, 2], [3, 4, 5]] -> memory 0 3 1 4 2 5
    const std::vector<int> col_major{0, 3, 1, 4, 2, 5};
    cnpy::save_npy(dir / "f.npy", col_major, {2, 3}, true);
    const auto a = cnpy::load_npy(dir / "f.npy");
    CHECK(a.fortran_order());
    CHECK(a.at<int>(0, 2) == 2);
    CHECK(a.at<int>(1, 0) == 3);
    CHECK(a.at<int>(1, 2) == 5);
}

TEST_CASE("typed access is checked") {
    temp_dir dir;
    cnpy::save_npy(dir / "d.npy", std::vector<double>{1, 2, 3});
    auto a = cnpy::load_npy(dir / "d.npy");
    CHECK_THROWS_AS((void)a.data<float>(), cnpy::dtype_error);
    CHECK_THROWS_AS((void)a.data<std::int64_t>(), cnpy::dtype_error);
    CHECK_THROWS_AS((void)a.at<double>(3), std::out_of_range);
    CHECK_THROWS_AS((void)a.at<double>(0, 0), std::out_of_range);
    a.mutable_data<double>()[0] = 42;
    CHECK(a.at<double>(0) == 42);
}

TEST_CASE("ndarray construction and copy semantics") {
    cnpy::ndarray z(cnpy::dtype::of<int>(), {2, 2});
    CHECK(z.to_vector<int>() == std::vector<int>{0, 0, 0, 0});
    CHECK(z.writable());

    auto a = cnpy::ndarray::from(std::vector<float>{1, 2, 3, 4}, {2, 2});
    const cnpy::ndarray shallow = a; // NOLINT: intentional
    cnpy::ndarray deep = a.copy();
    a.mutable_data<float>()[0] = 9;
    CHECK(shallow.at<float>(0, 0) == 9);
    CHECK(deep.at<float>(0, 0) == 1);
    CHECK(reinterpret_cast<std::uintptr_t>(deep.bytes().data()) % 64 == 0);

    CHECK_THROWS_AS((void)cnpy::ndarray::from(std::vector<float>{1, 2, 3}, {2, 2}), std::invalid_argument);
}

TEST_CASE("save_npy validates its input") {
    temp_dir dir;
    CHECK_THROWS_AS(cnpy::save_npy(dir / "x.npy", std::vector<int>{1, 2, 3}, {2, 2}), std::invalid_argument);
}

TEST_CASE("errors carry useful information") {
    temp_dir dir;
    SUBCASE("missing file") {
        try {
            (void)cnpy::load_npy(dir / "missing.npy");
            FAIL("expected io_error");
        } catch (const cnpy::io_error& e) {
            CHECK(e.path() == dir / "missing.npy");
            CHECK(e.code());
            CHECK(std::string(e.what()).find("missing.npy") != std::string::npos);
        }
    }
    SUBCASE("bad magic") {
        test::write_file(dir / "bad.npy", "not a numpy file at all");
        CHECK_THROWS_AS((void)cnpy::load_npy(dir / "bad.npy"), cnpy::format_error);
    }
    SUBCASE("truncated data") {
        cnpy::save_npy(dir / "t.npy", std::vector<double>(100));
        std::string f = test::read_file(dir / "t.npy");
        f.resize(f.size() - 8);
        test::write_file(dir / "t.npy", f);
        CHECK_THROWS_AS((void)cnpy::load_npy(dir / "t.npy"), cnpy::format_error);
    }
    SUBCASE("malformed header") {
        test::write_file(dir / "m.npy", test::handmade_npy("{'descr': '<f8', 'fortran_order': Maybe, 'shape': (1,), }",
                                                           std::string(8, '\0')));
        CHECK_THROWS_AS((void)cnpy::load_npy(dir / "m.npy"), cnpy::format_error);
    }
    SUBCASE("structured dtype") {
        test::write_file(dir / "r.npy",
                         test::handmade_npy("{'descr': [('a', '<i4'), ('b', '<f8')], 'fortran_order': False, "
                                            "'shape': (1,), }",
                                            std::string(12, '\0')));
        CHECK_THROWS_AS((void)cnpy::load_npy(dir / "r.npy"), cnpy::format_error);
    }
    SUBCASE("object dtype") {
        test::write_file(dir / "o.npy",
                         test::handmade_npy("{'descr': '|O', 'fortran_order': False, 'shape': (1,), }", "x"));
        CHECK_THROWS_AS((void)cnpy::load_npy(dir / "o.npy"), cnpy::format_error);
    }
}

TEST_CASE("headers written by other tools are accepted") {
    temp_dir dir;
    const auto v = values<std::int32_t>(6);
    SUBCASE("16-byte aligned header from cnpy 1.x, Python 2 long suffix, no trailing comma") {
        test::write_file(
            dir / "old.npy",
            test::handmade_npy("{'descr': '<i4', 'fortran_order': False, 'shape': (2L, 3L)}", test::raw_bytes(v)));
        const auto a = cnpy::load_npy(dir / "old.npy");
        CHECK(a.shape() == cnpy::shape_t{2, 3});
        CHECK(test::same(a.data<std::int32_t>(), v));
    }
    SUBCASE("keys in another order, double quotes") {
        test::write_file(dir / "q.npy", test::handmade_npy(R"({"shape": (6,), "fortran_order": False, "descr": "<i4"})",
                                                           test::raw_bytes(v)));
        CHECK(test::same(cnpy::load_npy(dir / "q.npy").data<std::int32_t>(), v));
    }
}

TEST_CASE("foreign byte order is converted on load") {
    temp_dir dir;
    const bool big_host = cnpy::native_byte_order == cnpy::byte_order::big;
    const auto v = values<std::int32_t>(5);
    const auto c = values<std::complex<double>>(4);
    test::write_file(dir / "be.npy", test::handmade_npy(std::string("{'descr': '") + (big_host ? "<" : ">") +
                                                            "i4', 'fortran_order': False, 'shape': (5,), }",
                                                        test::raw_bytes(v, true)));
    test::write_file(dir / "bec.npy", test::handmade_npy(std::string("{'descr': '") + (big_host ? "<" : ">") +
                                                             "c16', 'fortran_order': False, 'shape': (4,), }",
                                                         test::raw_bytes(c, true, sizeof(double))));
    const auto a = cnpy::load_npy(dir / "be.npy");
    CHECK(a.dtype().is_native());
    CHECK(test::same(a.data<std::int32_t>(), v));
    CHECK(test::same(cnpy::load_npy(dir / "bec.npy").data<std::complex<double>>(), c));
}

TEST_CASE("unicode file names") {
    temp_dir dir;
    const fs::path path = dir / fs::path(u8"größe_数据_😀.npy");
    cnpy::save_npy(path, std::vector<int>{7, 8});
    CHECK(fs::exists(path));
    CHECK(cnpy::load_npy(path).to_vector<int>() == std::vector<int>{7, 8});
}

TEST_CASE("version 2.0 header for very large dictionaries") {
    temp_dir dir;
    const cnpy::shape_t shape(25000, 1); // 25000 dimensions -> header > 64 KiB
    const double x = 3.0;
    cnpy::save_npy(dir / "v2.npy", std::span(&x, 1), shape);
    const std::string f = test::read_file(dir / "v2.npy");
    CHECK(f[6] == 2);
    const auto a = cnpy::load_npy(dir / "v2.npy");
    CHECK(a.ndim() == 25000);
    CHECK(a.data<double>()[0] == 3.0);
}

TEST_CASE("append_npy") {
    temp_dir dir;
    const auto path = dir / "app.npy";

    SUBCASE("grows axis 0 in place") {
        const auto v = values<double>(12);
        cnpy::append_npy(path, std::span(v).first(4), {1, 4}); // creates the file
        const auto size_before = test::read_file(path).size();
        cnpy::append_npy(path, std::span(v).subspan(4, 8), {2, 4});
        const auto a = cnpy::load_npy(path);
        CHECK(a.shape() == cnpy::shape_t{3, 4});
        CHECK(test::same(a.data<double>(), v));
        CHECK(test::read_file(path).size() == size_before + 8 * sizeof(double)); // header kept its size

        for (int i = 0; i < 200; ++i)
            cnpy::append_npy(path, std::span(v).first(4), {1, 4});
        CHECK(cnpy::load_npy(path).shape() == cnpy::shape_t{203, 4});
    }
    SUBCASE("1-d convenience") {
        cnpy::append_npy(path, std::vector<int>{1, 2});
        cnpy::append_npy(path, std::vector<int>{3});
        CHECK(cnpy::load_npy(path).to_vector<int>() == std::vector<int>{1, 2, 3});
    }
    SUBCASE("Fortran order grows the last axis") {
        cnpy::append_npy(path, std::vector<int>{0, 1}, {2, 1}, true);
        cnpy::append_npy(path, std::vector<int>{2, 3, 4, 5}, {2, 2}, true);
        const auto a = cnpy::load_npy(path);
        CHECK(a.shape() == cnpy::shape_t{2, 3});
        CHECK(a.at<int>(1, 2) == 5);
        CHECK(a.at<int>(0, 1) == 2);
    }
    SUBCASE("mismatches are rejected") {
        cnpy::save_npy(path, std::vector<double>(6), {2, 3});
        CHECK_THROWS_AS(cnpy::append_npy(path, std::vector<float>(3), {1, 3}), cnpy::dtype_error);
        CHECK_THROWS_AS(cnpy::append_npy(path, std::vector<double>(4), {1, 4}), cnpy::error);
        CHECK_THROWS_AS(cnpy::append_npy(path, std::vector<double>(3), {3}), cnpy::error);
        CHECK_THROWS_AS(cnpy::append_npy(path, std::vector<double>(3), {1, 3}, true), cnpy::error);
        CHECK(cnpy::load_npy(path).shape() == cnpy::shape_t{2, 3}); // untouched
    }
    SUBCASE("legacy header without growth room is rewritten") {
        const auto v = values<std::int16_t>(4);
        test::write_file(path, test::handmade_npy("{'descr': '<i2', 'fortran_order': False, 'shape': (2, 2), }",
                                                  test::raw_bytes(values<std::int16_t>(4)), 1));
        // Same number of digits: the header is rewritten in place.
        cnpy::append_npy(path, std::vector<std::int16_t>{1, 2}, {1, 2});
        CHECK(cnpy::load_npy(path).shape() == cnpy::shape_t{3, 2});
        // No spare room for more digits: goes through a temporary file.
        std::vector<std::int16_t> big(2 * 1000, 7);
        cnpy::append_npy(path, big, {1000, 2});
        const auto a = cnpy::load_npy(path);
        CHECK(a.shape() == cnpy::shape_t{1003, 2});
        CHECK(a.at<std::int16_t>(0, 1) == v[1]);
        CHECK(a.at<std::int16_t>(2, 1) == 2);
        CHECK(a.at<std::int16_t>(1002, 1) == 7);
        CHECK_FALSE(fs::exists(dir / "app.npy.cnpy-tmp"));
    }
}

#include "helpers.hpp"

#include <doctest/doctest.h>

using test::temp_dir;
using test::values;

TEST_CASE("map_npy read_only") {
    temp_dir dir;
    const auto v = values<double>(1000);
    cnpy::save_npy(dir / "m.npy", v, {10, 100});

    const auto a = cnpy::map_npy(dir / "m.npy");
    CHECK(a.is_mapped());
    CHECK_FALSE(a.writable());
    CHECK(a.shape() == cnpy::shape_t{10, 100});
    CHECK(test::same(a.data<double>(), v));
    CHECK(a.at<double>(9, 99) == v.back());
    // numpy-style headers are 64-byte multiples, so mapped data is 64-byte aligned
    CHECK(reinterpret_cast<std::uintptr_t>(a.bytes().data()) % 64 == 0);

    auto b = a; // shares the mapping
    CHECK_THROWS_AS((void)b.mutable_data<double>(), cnpy::error);
    CHECK_THROWS_AS((void)b.mutable_bytes(), cnpy::error);
}

TEST_CASE("map_npy outlives the mapped file handle and copies are independent") {
    temp_dir dir;
    cnpy::save_npy(dir / "m.npy", std::vector<int>{1, 2, 3});
    cnpy::ndarray owned;
    {
        const auto a = cnpy::map_npy(dir / "m.npy");
        owned = a.copy();
        CHECK_FALSE(owned.is_mapped());
        CHECK(owned.writable());
    }
    CHECK(owned.to_vector<int>() == std::vector<int>{1, 2, 3});
}

TEST_CASE("map_npy copy_on_write leaves the file untouched") {
    temp_dir dir;
    cnpy::save_npy(dir / "m.npy", std::vector<float>{1, 2, 3});
    {
        auto a = cnpy::map_npy(dir / "m.npy", cnpy::map_mode::copy_on_write);
        CHECK(a.writable());
        a.mutable_data<float>()[1] = 42;
        CHECK(a.at<float>(1) == 42);
    }
    CHECK(cnpy::load_npy(dir / "m.npy").to_vector<float>() == std::vector<float>{1, 2, 3});
}

TEST_CASE("map_npy read_write modifies the file") {
    temp_dir dir;
    cnpy::save_npy(dir / "m.npy", std::vector<std::int64_t>{1, 2, 3});
    {
        auto a = cnpy::map_npy(dir / "m.npy", cnpy::map_mode::read_write);
        a.mutable_data<std::int64_t>()[2] = -7;
    }
    CHECK(cnpy::load_npy(dir / "m.npy").to_vector<std::int64_t>() == std::vector<std::int64_t>{1, 2, -7});
}

TEST_CASE("map_npy keeps foreign byte order until as_native()") {
    temp_dir dir;
    const bool big_host = cnpy::native_byte_order == cnpy::byte_order::big;
    const auto v = values<std::uint32_t>(8);
    test::write_file(dir / "be.npy", test::handmade_npy(std::string("{'descr': '") + (big_host ? "<" : ">") +
                                                            "u4', 'fortran_order': False, 'shape': (8,), }",
                                                        test::raw_bytes(v, true), 64));
    const auto a = cnpy::map_npy(dir / "be.npy");
    CHECK_FALSE(a.dtype().is_native());
    CHECK_THROWS_AS((void)a.data<std::uint32_t>(), cnpy::dtype_error);
    const auto n = a.as_native();
    CHECK_FALSE(n.is_mapped());
    CHECK(test::same(n.data<std::uint32_t>(), v));
}

TEST_CASE("map_npy rejects truncated files") {
    temp_dir dir;
    cnpy::save_npy(dir / "t.npy", std::vector<double>(10));
    std::string f = test::read_file(dir / "t.npy");
    f.resize(f.size() - 1);
    test::write_file(dir / "t.npy", f);
    CHECK_THROWS_AS((void)cnpy::map_npy(dir / "t.npy"), cnpy::format_error);
}

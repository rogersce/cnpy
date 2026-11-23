#include "helpers.hpp"

#include <doctest/doctest.h>

#include <complex>
#include <cstdint>

using test::temp_dir;
using test::values;

namespace {
const cnpy::compression stored = cnpy::compression::none();
const cnpy::compression deflated = cnpy::compression::deflate();

/// Offset of the first byte of array data of the first member (stored, numpy-style header).
std::size_t first_member_data_offset(const std::string& zip) {
    const std::size_t name_len = test::le16(zip, 26);
    const std::size_t extra_len = test::le16(zip, 28);
    const std::size_t npy = 30 + name_len + extra_len;
    const std::size_t hlen = test::le16(zip, npy + 8);
    return npy + 10 + hlen;
}
} // namespace

TEST_CASE("npz_writer / npz_reader round-trip") {
    for (const auto c : {stored, deflated, cnpy::compression::deflate(1), cnpy::compression::deflate(9)}) {
        CAPTURE(c.level);
        temp_dir dir;
        const auto path = dir / "a.npz";
        const auto d = values<double>(60);
        const auto z = values<std::complex<float>>(7);
        const auto b = values<bool>(5);
        {
            cnpy::npz_writer w(path);
            w.add("d", d, {3, 4, 5}, c);
            w.add("z", z, c);
            w.add("b", b, c);
            w.add("empty", std::vector<int>{}, {0}, c);
            w.close();
        }
        const cnpy::npz_reader r(path);
        CHECK(r.size() == 4);
        CHECK(r.names() == std::vector<std::string>{"d", "z", "b", "empty"});
        CHECK(r.contains("z"));
        CHECK_FALSE(r.contains("q"));
        CHECK(r.load("d").shape() == cnpy::shape_t{3, 4, 5});
        CHECK(test::same(r.load("d").data<double>(), d));
        CHECK(test::same(r.load("z").data<std::complex<float>>(), z));
        CHECK(test::same(r.load("b").data<bool>(), b));
        CHECK(r.load("empty").empty());
        CHECK_THROWS_AS((void)r.load("q"), cnpy::error);

        const auto all = cnpy::load_npz(path);
        CHECK(all.size() == 4);
        CHECK(test::same(all.at("d").data<double>(), d));
        CHECK(test::same(cnpy::load_npz(path, "z").data<std::complex<float>>(), z));
    }
}

TEST_CASE("compression actually compresses") {
    temp_dir dir;
    const std::vector<std::int32_t> zeros(100000, 0);
    cnpy::save_npz(dir / "s.npz", "x", zeros);
    cnpy::save_npz(dir / "c.npz", "x", zeros, cnpy::npz_mode::truncate, deflated);
    CHECK(fs::file_size(dir / "c.npz") * 50 < fs::file_size(dir / "s.npz"));
    CHECK(cnpy::load_npz(dir / "c.npz", "x").to_vector<std::int32_t>() == zeros);
}

TEST_CASE("save_npz append mode and whole-map overload") {
    temp_dir dir;
    const auto path = dir / "m.npz";
    cnpy::save_npz(path, "a", std::vector<int>{1, 2, 3});                                // truncate
    cnpy::save_npz(path, "b", std::vector<float>{4, 5}, {2, 1}, cnpy::npz_mode::append); // stored
    cnpy::save_npz(path, "c", values<double>(1000), cnpy::npz_mode::append, deflated);
    {
        cnpy::npz_writer w(path, cnpy::npz_mode::append);
        CHECK_THROWS_AS(w.add("a", std::vector<int>{1}), std::invalid_argument); // duplicate
    }
    const auto all = cnpy::load_npz(path);
    REQUIRE(all.size() == 3);
    CHECK(all.at("a").to_vector<int>() == std::vector<int>{1, 2, 3});
    CHECK(all.at("b").shape() == cnpy::shape_t{2, 1});
    CHECK(test::same(all.at("c").data<double>(), values<double>(1000)));

    cnpy::save_npz(dir / "copy.npz", all, deflated);
    CHECK(cnpy::npz_reader(dir / "copy.npz").names() == std::vector<std::string>{"a", "b", "c"});
    CHECK(test::same(cnpy::load_npz(dir / "copy.npz", "c").data<double>(), values<double>(1000)));

    cnpy::save_npz(path, "only", std::vector<int>{1}); // truncate again
    CHECK(cnpy::npz_reader(path).size() == 1);
}

TEST_CASE("append to a missing archive creates it") {
    temp_dir dir;
    cnpy::save_npz(dir / "new.npz", "x", std::vector<int>{5}, cnpy::npz_mode::append);
    CHECK(cnpy::load_npz(dir / "new.npz", "x").to_vector<int>() == std::vector<int>{5});
}

TEST_CASE("npz_writer lifecycle") {
    temp_dir dir;
    SUBCASE("destructor finalises the archive") {
        {
            cnpy::npz_writer w(dir / "d.npz");
            w.add("x", std::vector<int>{1, 2});
        }
        CHECK(cnpy::load_npz(dir / "d.npz", "x").to_vector<int>() == std::vector<int>{1, 2});
    }
    SUBCASE("close is idempotent, add after close throws") {
        cnpy::npz_writer w(dir / "c.npz");
        w.close();
        w.close();
        CHECK_THROWS_AS(w.add("x", std::vector<int>{1}), cnpy::error);
        CHECK(cnpy::npz_reader(dir / "c.npz").size() == 0);
    }
    SUBCASE("move") {
        cnpy::npz_writer a(dir / "m.npz");
        a.add("x", std::vector<int>{1});
        cnpy::npz_writer b = std::move(a);
        b.add("y", std::vector<int>{2});
        b.close();
        CHECK(cnpy::npz_reader(dir / "m.npz").names() == std::vector<std::string>{"x", "y"});
    }
    SUBCASE("argument validation") {
        cnpy::npz_writer w(dir / "v.npz");
        CHECK_THROWS_AS(w.add("", std::vector<int>{1}), std::invalid_argument);
        CHECK_THROWS_AS(w.add("x", std::vector<int>{1}, cnpy::compression{10}), std::invalid_argument);
        CHECK_THROWS_AS(w.add("x", std::vector<int>{1, 2, 3}, {2, 2}), std::invalid_argument);
        w.add("x", std::vector<int>{1}); // the writer is still usable
        w.close();
        CHECK(cnpy::npz_reader(dir / "v.npz").names() == std::vector<std::string>{"x"});
    }
}

TEST_CASE("ndarray members keep dtype and memory order") {
    temp_dir dir;
    const auto f = cnpy::ndarray::from(std::vector<int>{0, 3, 1, 4, 2, 5}, {2, 3}, true);
    cnpy::npz_writer w(dir / "f.npz");
    w.add("f", f, deflated);
    w.close();
    const auto a = cnpy::load_npz(dir / "f.npz", "f");
    CHECK(a.fortran_order());
    CHECK(a.at<int>(1, 2) == 5);
}

TEST_CASE("UTF-8 member names") {
    temp_dir dir;
    const std::string name = "gr\xC3\xB6\xC3\x9F"
                             "e"; // "größe"
    cnpy::save_npz(dir / "u.npz", name, std::vector<int>{1});
    CHECK(cnpy::npz_reader(dir / "u.npz").names() == std::vector<std::string>{name});
}

TEST_CASE("CRC-32 detects corruption") {
    temp_dir dir;
    const auto path = dir / "crc.npz";
    cnpy::save_npz(path, "x", values<double>(100));
    std::string zip = test::read_file(path);
    const std::size_t at = first_member_data_offset(zip) + 17;
    zip[at] = static_cast<char>(zip[at] ^ 0x01);
    test::write_file(path, zip);

    CHECK_THROWS_AS((void)cnpy::load_npz(path, "x"), cnpy::format_error);
    cnpy::npz_read_options no_crc;
    no_crc.verify_crc = false;
    CHECK(cnpy::load_npz(path, "x", no_crc).size() == 100);
}

TEST_CASE("corrupt deflate stream is detected") {
    temp_dir dir;
    const auto path = dir / "bad.npz";
    cnpy::save_npz(path, "x", values<double>(5000), cnpy::npz_mode::truncate, deflated);
    std::string zip = test::read_file(path);
    for (std::size_t i = 200; i < 260; ++i)
        zip[i] = static_cast<char>(~zip[i]);
    test::write_file(path, zip);
    CHECK_THROWS_AS((void)cnpy::load_npz(path, "x"), cnpy::format_error);
}

TEST_CASE("not a zip file") {
    temp_dir dir;
    test::write_file(dir / "x.npz", std::string(100, 'x'));
    CHECK_THROWS_AS(cnpy::npz_reader(dir / "x.npz"), cnpy::format_error);
    test::write_file(dir / "tiny.npz", "PK");
    CHECK_THROWS_AS(cnpy::npz_reader(dir / "tiny.npz"), cnpy::format_error);
    CHECK_THROWS_AS(cnpy::npz_reader(dir / "missing.npz"), cnpy::io_error);
}

TEST_CASE("memory-mapped reader") {
    temp_dir dir;
    const auto path = dir / "m.npz";
    const auto bytes = values<std::uint8_t>(4096);
    const auto doubles = values<double>(333);
    {
        cnpy::npz_writer w(path);
        w.add("bytes", bytes);
        w.add("doubles", doubles);
        w.add("packed", doubles, deflated);
        w.close();
    }
    cnpy::npz_read_options opt;
    opt.memory_map = true;
    const cnpy::npz_reader r(path, opt);

    const auto a = r.load("bytes");
    CHECK(a.is_mapped()); // 1-byte elements are always aligned: zero-copy
    CHECK_FALSE(a.writable());
    CHECK(test::same(a.data<std::uint8_t>(), bytes));
    CHECK(test::same(r.load("doubles").data<double>(), doubles)); // mapped or copied, depending on alignment
    CHECK(test::same(r.load("packed").data<double>(), doubles));
    const auto all = cnpy::load_npz(path, opt);
    CHECK(all.size() == 3);
}

TEST_CASE("zip64: more than 65535 members") {
    temp_dir dir;
    const auto path = dir / "many.npz";
    constexpr std::size_t n = 65536 + 10;
    {
        cnpy::npz_writer w(path);
        for (std::size_t i = 0; i < n; ++i) {
            const auto v = static_cast<std::uint32_t>(i);
            w.add("a" + std::to_string(i), std::span(&v, 1), cnpy::shape_t{});
        }
        w.close();
    }
    const std::string zip = test::read_file(path);
    CHECK(zip.find(std::string("PK\x06\x06", 4)) != std::string::npos); // zip64 end record
    const cnpy::npz_reader r(path);
    CHECK(r.size() == n);
    CHECK(r.load("a0").at<std::uint32_t>() == 0);
    CHECK(r.load("a65540").at<std::uint32_t>() == 65540);

    cnpy::save_npz(path, "extra", std::vector<int>{1}, cnpy::npz_mode::append);
    CHECK(cnpy::npz_reader(path).size() == n + 1);
}

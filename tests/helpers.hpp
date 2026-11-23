#pragma once

#include <cnpy/cnpy.hpp>

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace test {

/// Unique scratch directory, removed on destruction.
class temp_dir {
public:
    temp_dir() {
        static std::atomic<int> counter{0};
        std::random_device rd;
        path_ = fs::temp_directory_path() /
                ("cnpy-test-" + std::to_string(rd()) + "-" + std::to_string(counter.fetch_add(1)));
        fs::create_directories(path_);
    }
    ~temp_dir() {
        std::error_code ec;
        fs::remove_all(path_, ec);
    }
    temp_dir(const temp_dir&) = delete;
    temp_dir& operator=(const temp_dir&) = delete;

    [[nodiscard]] fs::path operator/(const fs::path& name) const { return path_ / name; }
    [[nodiscard]] const fs::path& path() const noexcept { return path_; }

private:
    fs::path path_;
};

inline std::size_t le16(const std::string& s, std::size_t at) {
    return static_cast<std::size_t>(static_cast<unsigned char>(s[at])) |
           (static_cast<std::size_t>(static_cast<unsigned char>(s[at + 1])) << 8u);
}

inline std::string read_file(const fs::path& p) {
    std::ifstream in(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

inline void write_file(const fs::path& p, const std::string& content) {
    std::ofstream out(p, std::ios::binary | std::ios::trunc);
    out.write(content.data(), static_cast<std::streamsize>(content.size()));
}

/// std::vector<bool> is not contiguous, so tests use this for bool data.
class bool_vector {
public:
    bool_vector() = default;
    explicit bool_vector(std::size_t n) : p_(new bool[n]()), n_(n) {}
    bool_vector(const bool_vector& o) : bool_vector(o.n_) { std::copy(o.begin(), o.end(), begin()); }
    bool_vector& operator=(const bool_vector& o) {
        bool_vector tmp(o);
        std::swap(p_, tmp.p_);
        std::swap(n_, tmp.n_);
        return *this;
    }
    bool_vector(bool_vector&&) noexcept = default;
    bool_vector& operator=(bool_vector&&) noexcept = default;

    [[nodiscard]] bool* data() noexcept { return p_.get(); }
    [[nodiscard]] const bool* data() const noexcept { return p_.get(); }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] bool* begin() noexcept { return p_.get(); }
    [[nodiscard]] bool* end() noexcept { return p_.get() + n_; }
    [[nodiscard]] const bool* begin() const noexcept { return p_.get(); }
    [[nodiscard]] const bool* end() const noexcept { return p_.get() + n_; }
    bool& operator[](std::size_t i) noexcept { return p_[i]; }
    const bool& operator[](std::size_t i) const noexcept { return p_[i]; }

private:
    std::unique_ptr<bool[]> p_;
    std::size_t n_ = 0;
};

/// A contiguous container of T (std::vector, except for bool).
template <class T>
using vec = std::conditional_t<std::is_same_v<T, bool>, bool_vector, std::vector<T>>;

/// Deterministic test values, identical to tests/interop/interop.py.
template <cnpy::element T>
T value(std::size_t i) {
    if constexpr (std::is_same_v<T, bool>) {
        return i % 3 == 0;
    } else if constexpr (std::is_integral_v<T>) {
        const auto v = static_cast<long long>((i * 7 + 3) % 100);
        if constexpr (std::is_signed_v<T>)
            return static_cast<T>(v - 50);
        else
            return static_cast<T>(v);
    } else if constexpr (std::is_floating_point_v<T>) {
        return static_cast<T>(static_cast<double>(i) * 0.25 - 8.0);
    } else {
        using F = typename T::value_type;
        return T(static_cast<F>(static_cast<double>(i) * 0.25 - 8.0),
                 static_cast<F>(0.0 - static_cast<double>(i) * 0.5));
    }
}

template <cnpy::element T>
vec<T> values(std::size_t n) {
    vec<T> v(n);
    for (std::size_t i = 0; i < n; ++i)
        v[i] = value<T>(i);
    return v;
}

/// Element-wise equality of two ranges (works for spans, vectors and bool_vector).
template <class A, class B>
bool same(const A& a, const B& b) {
    return std::ranges::equal(a, b);
}

/// Builds a .npy file by hand, padded to `align` bytes like old writers (cnpy 1.x used 16).
inline std::string handmade_npy(const std::string& dict, const std::string& data, std::size_t align = 16) {
    std::string d = dict;
    const std::size_t total = 10 + d.size() + 1;
    d.append((align - total % align) % align, ' ');
    d.push_back('\n');
    std::string out = "\x93NUMPY";
    out.push_back('\x01');
    out.push_back('\x00');
    out.push_back(static_cast<char>(d.size() & 0xFFu));
    out.push_back(static_cast<char>(d.size() >> 8u));
    return out + d + data;
}

template <class T>
std::string raw_bytes(const std::vector<T>& v, bool swap = false, std::size_t unit = sizeof(T)) {
    std::string s(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
    if (swap)
        for (std::size_t i = 0; i + unit <= s.size(); i += unit)
            std::reverse(s.begin() + static_cast<std::ptrdiff_t>(i), s.begin() + static_cast<std::ptrdiff_t>(i + unit));
    return s;
}

} // namespace test

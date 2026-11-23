// Native file I/O and memory mapping.
#include "detail.hpp"

#include <algorithm>
#include <cerrno>
#include <limits>
#include <new>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace cnpy {

namespace {
std::string io_message(const std::string& what, const fs::path& path, std::error_code ec) {
    std::string msg = what + ": '" + detail::utf8(path) + "'";
    if (ec) msg += ": " + ec.message();
    return msg;
}
} // namespace

io_error::io_error(const std::string& what, fs::path path, std::error_code ec)
    : error(io_message(what, path, ec)), path_(std::move(path)), code_(ec) {}

namespace detail {

void throw_io(const std::string& what, const fs::path& path, std::error_code ec) {
    throw io_error(what, path, ec);
}

void throw_last_io_error(const std::string& what, const fs::path& path) {
#if defined(_WIN32)
    throw_io(what, path, std::error_code(static_cast<int>(::GetLastError()), std::system_category()));
#else
    throw_io(what, path, std::error_code(errno, std::generic_category()));
#endif
}

std::shared_ptr<std::byte> allocate(std::size_t nbytes) {
    constexpr std::align_val_t alignment{npy_align};
    auto* p = static_cast<std::byte*>(::operator new(std::max<std::size_t>(nbytes, 1), alignment));
    return {p, [](std::byte* q) { ::operator delete(q, alignment); }};
}

// Upper bound for a single read/write system call (Windows takes a DWORD, Linux caps at ~2 GiB).
constexpr std::size_t max_io_chunk = std::size_t{1} << 30;

#if defined(_WIN32)

namespace {
HANDLE as_handle(std::intptr_t h) noexcept {
    return reinterpret_cast<HANDLE>(h);
}
constexpr std::intptr_t invalid_handle = -1; // == INVALID_HANDLE_VALUE
} // namespace

file::file(const fs::path& path, open_mode mode) : path_(path) {
    DWORD access = GENERIC_READ;
    DWORD disposition = OPEN_EXISTING;
    switch (mode) {
    case open_mode::read: break;
    case open_mode::read_write: access |= GENERIC_WRITE; break;
    case open_mode::read_write_create:
        access |= GENERIC_WRITE;
        disposition = OPEN_ALWAYS;
        break;
    case open_mode::create:
        access |= GENERIC_WRITE;
        disposition = CREATE_ALWAYS;
        break;
    }
    HANDLE h = ::CreateFileW(path.c_str(), access, FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE, nullptr,
                             disposition, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (h == INVALID_HANDLE_VALUE) throw_last_io_error("cannot open file", path);
    handle_ = reinterpret_cast<std::intptr_t>(h);
}

void file::close() {
    if (handle_ != invalid_handle) {
        ::CloseHandle(as_handle(handle_));
        handle_ = invalid_handle;
    }
}

std::uint64_t file::size() const {
    LARGE_INTEGER sz;
    if (!::GetFileSizeEx(as_handle(handle_), &sz)) throw_last_io_error("cannot query file size", path_);
    return static_cast<std::uint64_t>(sz.QuadPart);
}

std::size_t file::read_at(std::uint64_t offset, std::span<std::byte> out) const {
    std::size_t done = 0;
    while (done < out.size()) {
        const auto chunk = static_cast<DWORD>(std::min(out.size() - done, max_io_chunk));
        OVERLAPPED ov{};
        const std::uint64_t pos = offset + done;
        ov.Offset = static_cast<DWORD>(pos & 0xFFFFFFFFu);
        ov.OffsetHigh = static_cast<DWORD>(pos >> 32);
        DWORD got = 0;
        if (!::ReadFile(as_handle(handle_), out.data() + done, chunk, &got, &ov)) {
            if (::GetLastError() == ERROR_HANDLE_EOF) break;
            throw_last_io_error("read failed", path_);
        }
        if (got == 0) break;
        done += got;
    }
    return done;
}

void file::write_at(std::uint64_t offset, std::span<const std::byte> data) {
    std::size_t done = 0;
    while (done < data.size()) {
        const auto chunk = static_cast<DWORD>(std::min(data.size() - done, max_io_chunk));
        OVERLAPPED ov{};
        const std::uint64_t pos = offset + done;
        ov.Offset = static_cast<DWORD>(pos & 0xFFFFFFFFu);
        ov.OffsetHigh = static_cast<DWORD>(pos >> 32);
        DWORD put = 0;
        if (!::WriteFile(as_handle(handle_), data.data() + done, chunk, &put, &ov) || put == 0)
            throw_last_io_error("write failed", path_);
        done += put;
    }
}

void file::truncate(std::uint64_t size) {
    FILE_END_OF_FILE_INFO info{};
    info.EndOfFile.QuadPart = static_cast<LONGLONG>(size);
    if (!::SetFileInformationByHandle(as_handle(handle_), FileEndOfFileInfo, &info, sizeof(info)))
        throw_last_io_error("cannot resize file", path_);
}

std::shared_ptr<std::byte> map_file(const file& f, map_mode mode, std::uint64_t size) {
    if (size == 0) throw format_error("cannot memory-map an empty file");
    DWORD protect = PAGE_READONLY;
    DWORD access = FILE_MAP_READ;
    if (mode == map_mode::copy_on_write) {
        protect = PAGE_WRITECOPY;
        access = FILE_MAP_COPY;
    } else if (mode == map_mode::read_write) {
        protect = PAGE_READWRITE;
        access = FILE_MAP_WRITE;
    }
    HANDLE mapping = ::CreateFileMappingW(as_handle(f.native_handle()), nullptr, protect, 0, 0, nullptr);
    if (mapping == nullptr) throw_last_io_error("cannot memory-map file", f.path());
    void* view = ::MapViewOfFile(mapping, access, 0, 0, static_cast<SIZE_T>(size));
    const DWORD err = ::GetLastError();
    ::CloseHandle(mapping); // the view keeps the mapping object alive
    if (view == nullptr)
        throw_io("cannot memory-map file", f.path(), std::error_code(static_cast<int>(err), std::system_category()));
    return {static_cast<std::byte*>(view), [](std::byte* p) { ::UnmapViewOfFile(p); }};
}

#else // POSIX

namespace {
constexpr std::intptr_t invalid_handle = -1;
int as_fd(std::intptr_t h) noexcept {
    return static_cast<int>(h);
}
} // namespace

file::file(const fs::path& path, open_mode mode) : path_(path) {
    int flags = O_CLOEXEC;
    switch (mode) {
    case open_mode::read: flags |= O_RDONLY; break;
    case open_mode::read_write: flags |= O_RDWR; break;
    case open_mode::read_write_create: flags |= O_RDWR | O_CREAT; break;
    case open_mode::create: flags |= O_RDWR | O_CREAT | O_TRUNC; break;
    }
    int fd = -1;
    do {
        fd = ::open(path.c_str(), flags, 0666);
    } while (fd < 0 && errno == EINTR);
    if (fd < 0) throw_last_io_error("cannot open file", path);
    handle_ = fd;
}

void file::close() {
    if (handle_ != invalid_handle) {
        ::close(as_fd(handle_));
        handle_ = invalid_handle;
    }
}

std::uint64_t file::size() const {
    struct stat st {};
    if (::fstat(as_fd(handle_), &st) != 0) throw_last_io_error("cannot query file size", path_);
    return static_cast<std::uint64_t>(st.st_size);
}

std::size_t file::read_at(std::uint64_t offset, std::span<std::byte> out) const {
    std::size_t done = 0;
    while (done < out.size()) {
        const std::size_t chunk = std::min(out.size() - done, max_io_chunk);
        const ssize_t got = ::pread(as_fd(handle_), out.data() + done, chunk, static_cast<off_t>(offset + done));
        if (got < 0) {
            if (errno == EINTR) continue;
            throw_last_io_error("read failed", path_);
        }
        if (got == 0) break;
        done += static_cast<std::size_t>(got);
    }
    return done;
}

void file::write_at(std::uint64_t offset, std::span<const std::byte> data) {
    std::size_t done = 0;
    while (done < data.size()) {
        const std::size_t chunk = std::min(data.size() - done, max_io_chunk);
        const ssize_t put = ::pwrite(as_fd(handle_), data.data() + done, chunk, static_cast<off_t>(offset + done));
        if (put < 0) {
            if (errno == EINTR) continue;
            throw_last_io_error("write failed", path_);
        }
        done += static_cast<std::size_t>(put);
    }
}

void file::truncate(std::uint64_t size) {
    int rc = 0;
    do {
        rc = ::ftruncate(as_fd(handle_), static_cast<off_t>(size));
    } while (rc != 0 && errno == EINTR);
    if (rc != 0) throw_last_io_error("cannot resize file", path_);
}

std::shared_ptr<std::byte> map_file(const file& f, map_mode mode, std::uint64_t size) {
    if (size == 0) throw format_error("cannot memory-map an empty file");
    if (size > std::numeric_limits<std::size_t>::max()) throw format_error("file too large to memory-map");
    int prot = PROT_READ;
    int flags = MAP_SHARED;
    if (mode == map_mode::copy_on_write) {
        prot |= PROT_WRITE;
        flags = MAP_PRIVATE;
    } else if (mode == map_mode::read_write) {
        prot |= PROT_WRITE;
    }
    const auto len = static_cast<std::size_t>(size);
    void* p = ::mmap(nullptr, len, prot, flags, as_fd(f.native_handle()), 0);
    if (p == MAP_FAILED) throw_last_io_error("cannot memory-map file", f.path());
    return {static_cast<std::byte*>(p), [len](std::byte* q) { ::munmap(q, len); }};
}

#endif

file::~file() {
    close();
}

file::file(file&& other) noexcept
    : path_(std::move(other.path_)), handle_(std::exchange(other.handle_, invalid_handle)) {}

file& file::operator=(file&& other) noexcept {
    if (this != &other) {
        close();
        path_ = std::move(other.path_);
        handle_ = std::exchange(other.handle_, invalid_handle);
    }
    return *this;
}

void file::read_exact_at(std::uint64_t offset, std::span<std::byte> out) const {
    if (read_at(offset, out) != out.size()) throw format_error("unexpected end of file in '" + utf8(path_) + "'");
}

} // namespace detail
} // namespace cnpy

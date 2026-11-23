// .npz (zip) reading and writing, with deflate and zip64 support.
#include "detail.hpp"

#include <zlib.h>

#include <algorithm>
#include <limits>
#include <optional>
#include <unordered_set>
#include <vector>

namespace cnpy {

namespace {

using detail::file;
using detail::header_info;
using detail::load_le;
using detail::open_mode;
using detail::store_le;

constexpr std::uint32_t sig_local = 0x04034b50;
constexpr std::uint32_t sig_central = 0x02014b50;
constexpr std::uint32_t sig_eocd = 0x06054b50;
constexpr std::uint32_t sig_eocd64 = 0x06064b50;
constexpr std::uint32_t sig_eocd64_locator = 0x07064b50;
constexpr std::uint16_t zip64_extra_id = 0x0001;
constexpr std::uint32_t u32_max = 0xFFFFFFFFu;
constexpr std::uint16_t u16_max = 0xFFFFu;
constexpr std::uint16_t method_stored = 0;
constexpr std::uint16_t method_deflate = 8;
constexpr std::uint16_t flag_encrypted = 0x0001;
constexpr std::uint16_t flag_utf8 = 0x0800;
constexpr std::uint16_t dos_date_1980_01_01 = (1 << 5) | 1;
constexpr std::size_t local_header_size = 30;
constexpr std::size_t central_header_size = 46;
constexpr std::size_t eocd_size = 22;
constexpr std::size_t stream_chunk = std::size_t{1} << 20;
constexpr std::size_t small_member = std::size_t{1} << 16;
constexpr std::size_t zlib_max_chunk = std::size_t{1} << 30; // zlib counts in 32-bit uInt

struct zip_entry {
    std::string name; ///< full member name, including ".npy"
    std::uint16_t flags = 0;
    std::uint16_t method = method_stored;
    std::uint16_t mod_time = 0;
    std::uint16_t mod_date = dos_date_1980_01_01;
    std::uint32_t crc = 0;
    std::uint64_t csize = 0;
    std::uint64_t usize = 0;
    std::uint64_t local_offset = 0;
};

struct central_directory {
    std::vector<zip_entry> entries;
    std::uint64_t offset = 0;
};

[[noreturn]] void bad_zip(const fs::path& p, const std::string& why) {
    throw format_error("'" + detail::utf8(p) + "' is not a valid .npz archive: " + why);
}

std::uint32_t crc32_of(std::uint32_t crc, std::span<const std::byte> data) noexcept {
    const auto* p = reinterpret_cast<const Bytef*>(data.data());
    std::size_t n = data.size();
    while (n > 0) { // crc32_z takes z_size_t, but stay portable to zlib builds where that is 32-bit
        const std::size_t chunk = std::min(n, zlib_max_chunk);
        crc = static_cast<std::uint32_t>(::crc32(crc, p, static_cast<uInt>(chunk)));
        p += chunk;
        n -= chunk;
    }
    return crc;
}

// ----------------------------------------------------------------------------- reading the directory

central_directory read_central_directory(const file& f, std::uint64_t file_size) {
    if (file_size < eocd_size) bad_zip(f.path(), "file too small");

    // The end-of-central-directory record sits within the last 22 + 65535 (max comment) bytes.
    const std::uint64_t tail_size = std::min<std::uint64_t>(file_size, eocd_size + u16_max);
    const std::uint64_t tail_start = file_size - tail_size;
    std::vector<std::byte> tail(static_cast<std::size_t>(tail_size));
    f.read_exact_at(tail_start, tail);

    std::optional<std::size_t> eocd;
    for (std::size_t i = tail.size() - eocd_size + 1; i-- > 0;) {
        if (load_le<std::uint32_t>(&tail[i]) == sig_eocd &&
            i + eocd_size + load_le<std::uint16_t>(&tail[i + 20]) <= tail.size()) {
            eocd = i;
            break;
        }
    }
    if (!eocd) bad_zip(f.path(), "end of central directory not found");

    const std::byte* e = &tail[*eocd];
    std::uint64_t n_entries = load_le<std::uint16_t>(e + 10);
    std::uint64_t cd_size = load_le<std::uint32_t>(e + 12);
    std::uint64_t cd_offset = load_le<std::uint32_t>(e + 16);
    if (load_le<std::uint16_t>(e + 4) != 0 || load_le<std::uint16_t>(e + 6) != 0)
        bad_zip(f.path(), "multi-disk archives are not supported");

    // zip64: a locator immediately precedes the classic record.
    const std::uint64_t eocd_abs = tail_start + *eocd;
    if (eocd_abs >= 20) {
        std::array<std::byte, 20> loc{};
        f.read_exact_at(eocd_abs - 20, loc);
        if (load_le<std::uint32_t>(loc.data()) == sig_eocd64_locator) {
            const auto eocd64_offset = load_le<std::uint64_t>(loc.data() + 8);
            std::array<std::byte, 56> rec{};
            if (eocd64_offset > file_size || file_size - eocd64_offset < rec.size())
                bad_zip(f.path(), "bad zip64 locator");
            f.read_exact_at(eocd64_offset, rec);
            if (load_le<std::uint32_t>(rec.data()) != sig_eocd64) bad_zip(f.path(), "bad zip64 end record");
            n_entries = load_le<std::uint64_t>(rec.data() + 32);
            cd_size = load_le<std::uint64_t>(rec.data() + 40);
            cd_offset = load_le<std::uint64_t>(rec.data() + 48);
        }
    }

    if (cd_offset > file_size || cd_size > file_size - cd_offset) bad_zip(f.path(), "central directory out of range");
    if (n_entries > cd_size / central_header_size) bad_zip(f.path(), "entry count exceeds directory size");

    std::vector<std::byte> cd(static_cast<std::size_t>(cd_size));
    f.read_exact_at(cd_offset, cd);

    central_directory dir;
    dir.offset = cd_offset;
    dir.entries.reserve(static_cast<std::size_t>(n_entries));
    std::size_t pos = 0;
    for (std::uint64_t k = 0; k < n_entries; ++k) {
        if (cd.size() - pos < central_header_size || load_le<std::uint32_t>(&cd[pos]) != sig_central)
            bad_zip(f.path(), "corrupt central directory");
        const std::byte* h = &cd[pos];
        zip_entry z;
        z.flags = load_le<std::uint16_t>(h + 8);
        z.method = load_le<std::uint16_t>(h + 10);
        z.mod_time = load_le<std::uint16_t>(h + 12);
        z.mod_date = load_le<std::uint16_t>(h + 14);
        z.crc = load_le<std::uint32_t>(h + 16);
        z.csize = load_le<std::uint32_t>(h + 20);
        z.usize = load_le<std::uint32_t>(h + 24);
        const std::size_t name_len = load_le<std::uint16_t>(h + 28);
        const std::size_t extra_len = load_le<std::uint16_t>(h + 30);
        const std::size_t comment_len = load_le<std::uint16_t>(h + 32);
        z.local_offset = load_le<std::uint32_t>(h + 42);
        if (cd.size() - pos - central_header_size < name_len + extra_len + comment_len)
            bad_zip(f.path(), "corrupt central directory");
        z.name.assign(reinterpret_cast<const char*>(h + central_header_size), name_len);

        // zip64 extended information: present fields replace the 0xFFFFFFFF placeholders, in order.
        const std::byte* x = h + central_header_size + name_len;
        const std::byte* x_end = x + extra_len;
        while (x_end - x >= 4) {
            const auto id = load_le<std::uint16_t>(x);
            const std::size_t len = load_le<std::uint16_t>(x + 2);
            x += 4;
            if (static_cast<std::size_t>(x_end - x) < len) break;
            if (id == zip64_extra_id) {
                const std::byte* v = x;
                auto take = [&](std::uint64_t& field) {
                    if (field != u32_max) return;
                    if (v + 8 > x + len) bad_zip(f.path(), "truncated zip64 field");
                    field = load_le<std::uint64_t>(v);
                    v += 8;
                };
                take(z.usize);
                take(z.csize);
                take(z.local_offset);
            }
            x += len;
        }
        pos += central_header_size + name_len + extra_len + comment_len;
        dir.entries.push_back(std::move(z));
    }
    return dir;
}

std::string_view strip_npy(std::string_view name) noexcept {
    if (name.size() >= 4 && name.substr(name.size() - 4) == ".npy") name.remove_suffix(4);
    return name;
}

std::size_t required_alignment(const dtype& t) noexcept {
    std::size_t a = t.itemsize;
    if (t.kind == 'c') a /= 2;
    if (t.kind == 'U') a = 4;
    if (t.kind == 'S' || t.kind == 'V' || t.kind == 'b' || a == 0) a = 1;
    a = std::min<std::size_t>(a, 16);
    return a & (~a + 1); // largest power of two dividing a
}

// ----------------------------------------------------------------------------- inflate

/// Streams a raw-deflate member out of a file (in 1 MiB chunks) or out of a memory mapping.
class inflater {
public:
    inflater(const file* f, std::uint64_t offset, std::uint64_t csize, std::span<const std::byte> mapped)
        : file_(f), offset_(offset), remaining_(csize), mapped_(mapped) {
        if (::inflateInit2(&zs_, -MAX_WBITS) != Z_OK) throw error("zlib inflateInit2 failed");
        if (file_ != nullptr) inbuf_.resize(static_cast<std::size_t>(std::min<std::uint64_t>(csize, stream_chunk)));
    }
    ~inflater() { ::inflateEnd(&zs_); }
    inflater(const inflater&) = delete;
    inflater& operator=(const inflater&) = delete;

    /// Fills `out` completely or throws format_error.
    void read(std::span<std::byte> out) {
        while (!out.empty()) {
            if (ended_) throw format_error("compressed .npz member is shorter than its header claims");
            const std::size_t chunk = std::min(out.size(), zlib_max_chunk);
            const std::size_t produced = step(out.first(chunk));
            crc_ = crc32_of(crc_, out.first(produced));
            total_ += produced;
            out = out.subspan(produced);
        }
    }

    /// Checks that the deflate stream ends exactly here.
    void finish() {
        std::byte dummy{};
        while (!ended_) {
            if (step({&dummy, 1}) != 0) throw format_error("compressed .npz member has trailing data");
        }
    }

    [[nodiscard]] std::uint32_t crc() const noexcept { return crc_; }
    [[nodiscard]] std::uint64_t total_out() const noexcept { return total_; }

private:
    std::size_t step(std::span<std::byte> out) {
        if (zs_.avail_in == 0) refill();
        zs_.next_out = reinterpret_cast<Bytef*>(out.data());
        zs_.avail_out = static_cast<uInt>(out.size());
        const int rc = ::inflate(&zs_, Z_NO_FLUSH);
        const std::size_t produced = out.size() - zs_.avail_out;
        if (rc == Z_STREAM_END) {
            ended_ = true;
        } else if (rc == Z_BUF_ERROR) {
            if (zs_.avail_in == 0 && remaining_ == 0) throw format_error("compressed .npz member is truncated");
        } else if (rc != Z_OK) {
            throw format_error(std::string("corrupt compressed .npz member: ") + (zs_.msg ? zs_.msg : "zlib error"));
        }
        return produced;
    }

    void refill() {
        if (remaining_ == 0) return;
        const std::size_t n = static_cast<std::size_t>(
            std::min<std::uint64_t>(remaining_, file_ != nullptr ? inbuf_.size() : zlib_max_chunk));
        const std::byte* src = nullptr;
        if (file_ != nullptr) {
            file_->read_exact_at(offset_, std::span(inbuf_).first(n));
            src = inbuf_.data();
        } else {
            src = mapped_.data() + (mapped_.size() - remaining_);
        }
        offset_ += n;
        remaining_ -= n;
        zs_.next_in = reinterpret_cast<z_const Bytef*>(const_cast<std::byte*>(src));
        zs_.avail_in = static_cast<uInt>(n);
    }

    z_stream zs_{};
    const file* file_;
    std::uint64_t offset_;
    std::uint64_t remaining_;
    std::span<const std::byte> mapped_;
    std::vector<std::byte> inbuf_;
    std::uint32_t crc_ = 0;
    std::uint64_t total_ = 0;
    bool ended_ = false;
};

// ----------------------------------------------------------------------------- writing records

bool needs_zip64(const zip_entry& z) noexcept {
    return z.usize >= u32_max || z.csize >= u32_max || z.local_offset >= u32_max;
}

std::string local_header(const zip_entry& z, bool zip64) {
    std::string h;
    h.reserve(local_header_size + z.name.size() + 20);
    store_le(h, sig_local);
    store_le(h, static_cast<std::uint16_t>(zip64 ? 45 : 20));
    store_le(h, z.flags);
    store_le(h, z.method);
    store_le(h, z.mod_time);
    store_le(h, z.mod_date);
    store_le(h, z.crc);
    store_le(h, zip64 ? u32_max : static_cast<std::uint32_t>(z.csize));
    store_le(h, zip64 ? u32_max : static_cast<std::uint32_t>(z.usize));
    store_le(h, static_cast<std::uint16_t>(z.name.size()));
    store_le(h, static_cast<std::uint16_t>(zip64 ? 20 : 0));
    h += z.name;
    if (zip64) {
        store_le(h, zip64_extra_id);
        store_le(h, std::uint16_t{16});
        store_le(h, z.usize);
        store_le(h, z.csize);
    }
    return h;
}

void append_central_header(std::string& out, const zip_entry& z) {
    std::string extra;
    if (z.usize >= u32_max) store_le(extra, z.usize);
    if (z.csize >= u32_max) store_le(extra, z.csize);
    if (z.local_offset >= u32_max) store_le(extra, z.local_offset);
    if (!extra.empty()) {
        std::string field;
        store_le(field, zip64_extra_id);
        store_le(field, static_cast<std::uint16_t>(extra.size()));
        extra = field + extra;
    }
    const bool zip64 = needs_zip64(z);
    store_le(out, sig_central);
    store_le(out, std::uint16_t{45}); // version made by: 4.5, MS-DOS attribute compatibility
    store_le(out, static_cast<std::uint16_t>(zip64 ? 45 : 20));
    store_le(out, z.flags);
    store_le(out, z.method);
    store_le(out, z.mod_time);
    store_le(out, z.mod_date);
    store_le(out, z.crc);
    store_le(out, static_cast<std::uint32_t>(std::min<std::uint64_t>(z.csize, u32_max)));
    store_le(out, static_cast<std::uint32_t>(std::min<std::uint64_t>(z.usize, u32_max)));
    store_le(out, static_cast<std::uint16_t>(z.name.size()));
    store_le(out, static_cast<std::uint16_t>(extra.size()));
    store_le(out, std::uint16_t{0}); // comment length
    store_le(out, std::uint16_t{0}); // disk number
    store_le(out, std::uint16_t{0}); // internal attributes
    store_le(out, std::uint32_t{0}); // external attributes
    store_le(out, static_cast<std::uint32_t>(std::min<std::uint64_t>(z.local_offset, u32_max)));
    out += z.name;
    out += extra;
}

bool is_ascii(std::string_view s) noexcept {
    return std::all_of(s.begin(), s.end(), [](char c) { return static_cast<unsigned char>(c) < 0x80; });
}

} // namespace

// ============================================================================= npz_reader

struct npz_reader::impl {
    fs::path path;
    npz_read_options options;
    std::optional<file> f;           // used when not memory-mapped
    std::shared_ptr<std::byte> base; // whole-file mapping
    std::uint64_t file_size = 0;
    std::vector<zip_entry> entries;
    std::map<std::string, std::size_t, std::less<>> index; // stripped name -> entry (last one wins)

    [[nodiscard]] std::span<const std::byte> mapped(std::uint64_t offset, std::uint64_t n) const {
        if (offset > file_size || n > file_size - offset) bad_zip(path, "member data out of range");
        return {base.get() + offset, static_cast<std::size_t>(n)};
    }

    void read(std::uint64_t offset, std::span<std::byte> out) const {
        if (base) {
            const auto src = mapped(offset, out.size());
            std::copy(src.begin(), src.end(), out.begin());
        } else {
            f->read_exact_at(offset, out);
        }
    }

    [[nodiscard]] std::uint64_t data_offset(const zip_entry& z) const {
        std::array<std::byte, local_header_size> lh{};
        if (z.local_offset > file_size || file_size - z.local_offset < lh.size()) bad_zip(path, "bad member offset");
        read(z.local_offset, lh);
        if (load_le<std::uint32_t>(lh.data()) != sig_local) bad_zip(path, "bad local file header");
        const std::uint64_t off = z.local_offset + local_header_size + load_le<std::uint16_t>(lh.data() + 26) +
                                  load_le<std::uint16_t>(lh.data() + 28);
        if (off > file_size || z.csize > file_size - off) bad_zip(path, "member data out of range");
        return off;
    }

    void verify(const zip_entry& z, std::uint32_t crc) const {
        if (crc != z.crc) throw format_error("CRC-32 mismatch for '" + z.name + "' in '" + detail::utf8(path) + "'");
    }

    [[nodiscard]] static ndarray finish(const header_info& h, std::shared_ptr<std::byte> buffer) {
        dtype type = h.type;
        if (!type.is_native()) {
            detail::byteswap({buffer.get(), h.nbytes}, type);
            type = type.native();
        }
        return detail::array_access::make(type, h.shape, h.fortran_order, std::move(buffer), true, false);
    }

    [[nodiscard]] ndarray load_stored(const zip_entry& z, std::uint64_t off) const {
        if (z.csize != z.usize) bad_zip(path, "stored member '" + z.name + "' has inconsistent sizes");
        header_info h;
        if (base) {
            h = detail::read_header(mapped(off, z.usize));
        } else {
            h = detail::read_header(*f, off, z.usize);
        }
        if (z.usize != h.header_size + h.nbytes)
            throw format_error("member '" + z.name + "' size does not match its .npy header");

        const std::uint64_t data_off = off + h.header_size;
        if (base) {
            std::byte* p = base.get() + data_off;
            if (reinterpret_cast<std::uintptr_t>(p) % required_alignment(h.type) == 0) {
                // Zero-copy view into the mapping (read-only, CRC not checked to keep it lazy).
                return detail::array_access::make(h.type, h.shape, h.fortran_order, std::shared_ptr<std::byte>(base, p),
                                                  false, true);
            }
        }
        auto buffer = detail::allocate(h.nbytes);
        read(data_off, {buffer.get(), h.nbytes});
        if (options.verify_crc) {
            std::vector<std::byte> hb(h.header_size);
            read(off, hb);
            verify(z, crc32_of(crc32_of(0, hb), {buffer.get(), h.nbytes}));
        }
        return finish(h, std::move(buffer));
    }

    [[nodiscard]] ndarray load_deflated(const zip_entry& z, std::uint64_t off) const {
        inflater inf(base ? nullptr : &*f, off, z.csize, base ? mapped(off, z.csize) : std::span<const std::byte>{});

        std::array<std::byte, detail::npy_prefix_v2> first{};
        const auto first_n = static_cast<std::size_t>(std::min<std::uint64_t>(z.usize, first.size()));
        inf.read(std::span(first).first(first_n));
        const detail::header_prefix pre = detail::parse_header_prefix(std::span(first).first(first_n));
        if (pre.header_size > z.usize) bad_zip(path, "truncated .npy header in '" + z.name + "'");

        std::vector<std::byte> header(pre.header_size);
        std::copy_n(first.begin(), first_n, header.begin());
        inf.read(std::span(header).subspan(first_n));
        const header_info h = detail::read_header(header);
        if (z.usize != h.header_size + h.nbytes)
            throw format_error("member '" + z.name + "' size does not match its .npy header");

        auto buffer = detail::allocate(h.nbytes);
        inf.read({buffer.get(), h.nbytes});
        inf.finish();
        if (options.verify_crc) verify(z, inf.crc());
        return finish(h, std::move(buffer));
    }

    [[nodiscard]] ndarray load(const zip_entry& z) const {
        if ((z.flags & flag_encrypted) != 0) throw format_error("encrypted .npz members are not supported");
        const std::uint64_t off = data_offset(z);
        switch (z.method) {
        case method_stored: return load_stored(z, off);
        case method_deflate: return load_deflated(z, off);
        default: throw format_error("unsupported zip compression method " + std::to_string(z.method));
        }
    }
};

npz_reader::npz_reader(const fs::path& path, npz_read_options options) : impl_(std::make_unique<impl>()) {
    impl_->path = path;
    impl_->options = options;
    file f(path, open_mode::read);
    impl_->file_size = f.size();
    central_directory dir = read_central_directory(f, impl_->file_size);
    impl_->entries = std::move(dir.entries);
    for (std::size_t i = 0; i < impl_->entries.size(); ++i)
        impl_->index.insert_or_assign(std::string(strip_npy(impl_->entries[i].name)), i);
    if (options.memory_map) {
        impl_->base = detail::map_file(f, map_mode::read_only, impl_->file_size);
    } else {
        impl_->f.emplace(std::move(f));
    }
}

npz_reader::~npz_reader() = default;
npz_reader::npz_reader(npz_reader&&) noexcept = default;
npz_reader& npz_reader::operator=(npz_reader&&) noexcept = default;

std::vector<std::string> npz_reader::names() const {
    std::vector<std::string> out;
    out.reserve(impl_->index.size());
    for (std::size_t i = 0; i < impl_->entries.size(); ++i) {
        const auto name = strip_npy(impl_->entries[i].name);
        if (const auto it = impl_->index.find(name); it != impl_->index.end() && it->second == i)
            out.emplace_back(name);
    }
    return out;
}

bool npz_reader::contains(std::string_view name) const {
    return impl_->index.contains(name);
}

std::size_t npz_reader::size() const noexcept {
    return impl_->index.size();
}

ndarray npz_reader::load(std::string_view name) const {
    const auto it = impl_->index.find(name);
    if (it == impl_->index.end())
        throw error("no array named '" + std::string(name) + "' in '" + detail::utf8(impl_->path) + "'");
    return impl_->load(impl_->entries[it->second]);
}

npz_t npz_reader::load_all() const {
    npz_t out;
    for (const auto& [name, i] : impl_->index)
        out.emplace(name, impl_->load(impl_->entries[i]));
    return out;
}

npz_t load_npz(const fs::path& path, npz_read_options options) {
    return npz_reader(path, options).load_all();
}

ndarray load_npz(const fs::path& path, std::string_view name, npz_read_options options) {
    return npz_reader(path, options).load(name);
}

// ============================================================================= npz_writer

struct npz_writer::impl {
    file f;
    std::vector<zip_entry> entries;
    std::unordered_set<std::string> names; // for O(1) duplicate detection
    std::uint64_t pos = 0;
    bool closed = false;

    explicit impl(file&& file_) : f(std::move(file_)) {}

    void write(std::span<const std::byte> data) {
        f.write_at(pos, data);
        pos += data.size();
    }
    void write(std::string_view s) { write(detail::as_bytes(s)); }

    void add_stored(zip_entry& z, std::string_view header, std::span<const std::byte> bytes) {
        z.crc = crc32_of(crc32_of(0, detail::as_bytes(header)), bytes);
        z.csize = z.usize;
        std::string head = local_header(z, z.usize >= u32_max);
        head += header;
        if (bytes.size() <= small_member) { // one system call for small arrays
            head.append(reinterpret_cast<const char*>(bytes.data()), bytes.size());
            write(head);
        } else {
            write(head);
            write(bytes);
        }
    }

    void add_deflated(zip_entry& z, std::string_view header, std::span<const std::byte> bytes, int level) {
        // Decide on zip64 before knowing the compressed size; deflate can expand data very slightly.
        const bool zip64 = z.usize >= u32_max - (z.usize / 1000 + (1u << 16));
        const std::string placeholder = local_header(z, zip64);
        write(placeholder);
        const std::uint64_t data_start = pos;

        z_stream zs{};
        if (::deflateInit2(&zs, level, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
            throw error("zlib deflateInit2 failed");
        struct guard {
            z_stream* s;
            ~guard() { ::deflateEnd(s); }
        } g{&zs};

        std::vector<std::byte> out(stream_chunk);
        auto pump = [&](std::span<const std::byte> in, bool last) {
            do {
                const std::size_t take = std::min(in.size(), zlib_max_chunk);
                zs.next_in = reinterpret_cast<z_const Bytef*>(const_cast<std::byte*>(in.data()));
                zs.avail_in = static_cast<uInt>(take);
                const bool final_chunk = last && take == in.size();
                int rc = Z_OK;
                do {
                    zs.next_out = reinterpret_cast<Bytef*>(out.data());
                    zs.avail_out = static_cast<uInt>(out.size());
                    rc = ::deflate(&zs, final_chunk ? Z_FINISH : Z_NO_FLUSH);
                    if (rc == Z_STREAM_ERROR) throw error("zlib deflate failed");
                    write(std::span(out).first(out.size() - zs.avail_out));
                } while (zs.avail_out == 0 || (final_chunk && rc != Z_STREAM_END));
                in = in.subspan(take);
            } while (!in.empty());
        };
        z.crc = crc32_of(crc32_of(0, detail::as_bytes(header)), bytes);
        pump(detail::as_bytes(header), false);
        pump(bytes, true);

        z.csize = pos - data_start;
        if (!zip64 && z.csize >= u32_max) throw error("internal error: deflate output exceeded zip32 limits");
        f.write_at(z.local_offset, detail::as_bytes(local_header(z, zip64)));
    }

    void close() {
        if (closed) return;
        const std::uint64_t cd_offset = pos;
        std::string cd;
        for (const auto& z : entries)
            append_central_header(cd, z);
        write(cd);
        const std::uint64_t cd_size = cd.size();
        const std::uint64_t n = entries.size();

        std::string tail;
        if (n >= u16_max || cd_size >= u32_max || cd_offset >= u32_max) {
            const std::uint64_t eocd64_offset = pos;
            store_le(tail, sig_eocd64);
            store_le(tail, std::uint64_t{44}); // size of the remaining record
            store_le(tail, std::uint16_t{45});
            store_le(tail, std::uint16_t{45});
            store_le(tail, std::uint32_t{0});
            store_le(tail, std::uint32_t{0});
            store_le(tail, n);
            store_le(tail, n);
            store_le(tail, cd_size);
            store_le(tail, cd_offset);
            store_le(tail, sig_eocd64_locator);
            store_le(tail, std::uint32_t{0});
            store_le(tail, eocd64_offset);
            store_le(tail, std::uint32_t{1});
        }
        store_le(tail, sig_eocd);
        store_le(tail, std::uint16_t{0});
        store_le(tail, std::uint16_t{0});
        store_le(tail, static_cast<std::uint16_t>(std::min<std::uint64_t>(n, u16_max)));
        store_le(tail, static_cast<std::uint16_t>(std::min<std::uint64_t>(n, u16_max)));
        store_le(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(cd_size, u32_max)));
        store_le(tail, static_cast<std::uint32_t>(std::min<std::uint64_t>(cd_offset, u32_max)));
        store_le(tail, std::uint16_t{0});
        write(tail);
        f.truncate(pos); // appending may leave stale bytes of the old directory behind
        f.close();
        closed = true;
    }
};

npz_writer::npz_writer(const fs::path& path, npz_mode mode) {
    if (mode == npz_mode::append) {
        file f(path, open_mode::read_write_create);
        const std::uint64_t size = f.size();
        impl_ = std::make_unique<impl>(std::move(f));
        if (size != 0) {
            central_directory dir = read_central_directory(impl_->f, size);
            impl_->entries = std::move(dir.entries);
            for (const auto& e : impl_->entries)
                impl_->names.insert(e.name);
            impl_->pos = dir.offset; // new members overwrite the old central directory
        }
    } else {
        impl_ = std::make_unique<impl>(file(path, open_mode::create));
    }
}

npz_writer::~npz_writer() {
    if (impl_ && !impl_->closed) {
        try {
            impl_->close();
        } catch (...) { // destructors must not throw; call close() to observe errors
        }
    }
}

npz_writer::npz_writer(npz_writer&&) noexcept = default;

npz_writer& npz_writer::operator=(npz_writer&& other) noexcept {
    if (this != &other) {
        if (impl_ && !impl_->closed) {
            try {
                impl_->close();
            } catch (...) {
            }
        }
        impl_ = std::move(other.impl_);
    }
    return *this;
}

void npz_writer::add(std::string_view name, std::span<const std::byte> bytes, const dtype& type,
                     std::span<const std::size_t> shape, bool fortran_order, compression c) {
    if (!impl_ || impl_->closed) throw error("npz_writer: archive is already closed");
    if (c.level < 0 || c.level > 9) throw std::invalid_argument("npz_writer: compression level must be 0..9");
    if (name.empty()) throw std::invalid_argument("npz_writer: array name must not be empty");
    detail::check_bytes(bytes, type, shape);

    zip_entry z;
    z.name = std::string(name) + ".npy";
    if (z.name.size() > u16_max) throw std::invalid_argument("npz_writer: array name too long");
    if (impl_->names.contains(z.name))
        throw std::invalid_argument("npz_writer: duplicate array name '" + std::string(name) + "'");
    if (!is_ascii(z.name)) z.flags |= flag_utf8;
    z.method = c.level == 0 ? method_stored : method_deflate;
    z.local_offset = impl_->pos;

    const std::string header = detail::make_header(type, shape, fortran_order);
    z.usize = header.size() + bytes.size();

    try {
        if (z.method == method_stored) {
            impl_->add_stored(z, header, bytes);
        } else {
            impl_->add_deflated(z, header, bytes, c.level);
        }
    } catch (...) {
        impl_->pos = z.local_offset; // drop the partial member; the directory is written over it
        throw;
    }
    impl_->names.insert(z.name);
    impl_->entries.push_back(std::move(z));
}

void npz_writer::add(std::string_view name, const ndarray& array, compression c) {
    add(name, array.bytes(), array.dtype(), array.shape(), array.fortran_order(), c);
}

void npz_writer::close() {
    if (impl_) impl_->close();
}

void save_npz(const fs::path& path, const npz_t& arrays, compression c) {
    npz_writer w(path, npz_mode::truncate);
    for (const auto& [name, array] : arrays)
        w.add(name, array, c);
    w.close();
}

} // namespace cnpy

#include "kmer_index.hpp"
#include "../io/fasta.hpp"
#include "../util/sequence.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// ---------------------------------------------------------------------------
// File magic
// ---------------------------------------------------------------------------
static constexpr char MAGIC[8] = {'T','E','B','_','I','D','X','1'};

// ---------------------------------------------------------------------------
// Base encoding table: A=0, C=1, G=2, T=3; N/unknown = 0xFF (invalid).
// Matches GenomeStorage encoding so packed data is byte-for-byte compatible.
// ---------------------------------------------------------------------------
static constexpr auto make_base_enc() noexcept {
    struct T { uint8_t v[256]; };
    T t{};
    for (auto& x : t.v) x = 0xFFu;
    t.v[static_cast<uint8_t>('A')] = t.v[static_cast<uint8_t>('a')] = 0;
    t.v[static_cast<uint8_t>('C')] = t.v[static_cast<uint8_t>('c')] = 1;
    t.v[static_cast<uint8_t>('G')] = t.v[static_cast<uint8_t>('g')] = 2;
    t.v[static_cast<uint8_t>('T')] = t.v[static_cast<uint8_t>('t')] = 3;
    return t;
}
static constexpr auto BASE_ENC = make_base_enc();

// ---------------------------------------------------------------------------
// Little-endian helpers (safe with any pointer alignment).
// ---------------------------------------------------------------------------
static inline uint32_t load_le32(const uint8_t* p) noexcept {
    uint32_t v; std::memcpy(&v, p, 4); return v;
}
static inline uint64_t load_le64(const uint8_t* p) noexcept {
    uint64_t v; std::memcpy(&v, p, 8); return v;
}

// ---------------------------------------------------------------------------
// KmerIndex constructor / destructor
// ---------------------------------------------------------------------------
KmerIndex::KmerIndex(size_t k) : k_(k), table_() {}

KmerIndex::~KmerIndex() {
    if (mmap_ptr_ && mmap_ptr_ != MAP_FAILED) {
        munmap(mmap_ptr_, mmap_size_);
        mmap_ptr_ = nullptr;
    }
    if (mmap_fd_ >= 0) {
        close(mmap_fd_);
        mmap_fd_ = -1;
    }
}

// ---------------------------------------------------------------------------
// build()  —  two-pass CSR build
//
// Binary k-mer encoding mirrors GenomeStorage: A=0, C=1, G=2, T=3.
// The running key is maintained as:
//   key = (oldest_base << 2*(k-1)) | ... | (newest_base << 0)
// Rolling step (shift in new base at LSB, old base naturally falls off the
// top after applying kmer_mask):
//   new_key = ((key << 2) | new_base_bits) & kmer_mask
// ---------------------------------------------------------------------------
// max_freq_ filter discards highly repetitive k-mers
void KmerIndex::build(const std::string& fasta_path) {
    // Parse → encode → free the raw string before any table work.
    {
        std::string seq = read_fasta_sequence(fasta_path);
        genome_.encode(std::move(seq));
    }

    packed_ptr_  = genome_.packed_data();
    nmask_ptr_   = genome_.n_mask_data();
    genome_size_ = genome_.size();

    const uint32_t len = genome_size_;
    if (len < static_cast<uint32_t>(k_)) {
        std::cerr << "[kmer_index] genome (" << len
                  << " bp) shorter than k=" << k_ << " — skipping\n";
        return;
    }

    const uint64_t kmer_mask = (k_ < 32)
        ? ((uint64_t{1} << (2 * k_)) - 1)
        : ~uint64_t{0};
    const uint32_t last_start = len - static_cast<uint32_t>(k_);
    constexpr uint32_t REPORT_STEP = 100'000'000u;

    const uint8_t*  pk = packed_ptr_;
    const uint64_t* nm = nmask_ptr_;
    auto packed_base = [pk](uint32_t pos) noexcept -> uint8_t {
        return (pk[pos / 4] >> ((pos % 4) * 2)) & 3u;
    };
    auto is_n_fast = [nm](uint32_t pos) noexcept -> bool {
        return (nm[pos / 64] >> (pos % 64)) & uint64_t{1};
    };

    // =========================================================================
    // Phase 1: count occurrences of each k-mer.
    // =========================================================================
    std::cerr << "[kmer_index] Phase 1/2: counting k-mers...\n";

    ankerl::unordered_dense::map<uint64_t, uint32_t> count_map;
    {
        uint32_t n_count = 0;
        for (uint32_t i = 0; i < static_cast<uint32_t>(k_); ++i)
            if (is_n_fast(i)) ++n_count;

        uint64_t key = 0;
        for (uint32_t i = 0; i < static_cast<uint32_t>(k_); ++i)
            key = ((key << 2) | packed_base(i)) & kmer_mask;

        if (n_count == 0) ++count_map[key];

        uint32_t next_report = REPORT_STEP;
        for (uint32_t i = 1; i <= last_start; ++i) {
            if (is_n_fast(i - 1))                             --n_count;
            if (is_n_fast(i + static_cast<uint32_t>(k_) - 1)) ++n_count;
            key = ((key << 2) | packed_base(i + static_cast<uint32_t>(k_) - 1))
                  & kmer_mask;
            if (n_count == 0) ++count_map[key];

            if (i >= next_report) {
                std::cerr << "[kmer_index] phase 1: " << (i / 1'000'000u)
                          << " / " << (len / 1'000'000u)
                          << " Mbp  distinct=" << count_map.size() << "\n";
                next_report += REPORT_STEP;
            }
        }
    }
    std::cerr << "[kmer_index] Phase 1 done. distinct_kmers="
              << count_map.size() << "\n";

    // =========================================================================
    // Transition: compute offsets, build table_, free count_map, alloc flat_.
    // =========================================================================
    uint64_t total_positions = 0;
    uint64_t filtered_kmers  = 0;
    for (const auto& [k, c] : count_map) {
        if (max_freq_ > 0 && c > max_freq_) { ++filtered_kmers; continue; }
        total_positions += c;
    }
    if (filtered_kmers > 0)
        std::cerr << "[kmer_index] Filtered " << filtered_kmers
                  << " k-mers with >" << max_freq_ << " occurrences\n";

    std::cerr << "[kmer_index] Building CSR structure ("
              << (total_positions * sizeof(uint32_t) / (1024ULL * 1024)) << " MB flat)...\n";

    table_.clear();
    table_.reserve(count_map.size() - filtered_kmers);
    {
        uint64_t offset = 0;
        for (const auto& [k, c] : count_map) {
            if (max_freq_ > 0 && c > max_freq_) continue;
            // second starts at 0 and acts as fill cursor in phase 2.
            table_.emplace(k, std::make_pair(static_cast<uint32_t>(offset),
                                             uint32_t{0}));
            offset += c;
        }
    }

    // Free count_map before allocating flat_ (its memory is now free).
    { ankerl::unordered_dense::map<uint64_t, uint32_t> tmp; count_map.swap(tmp); }

    flat_.clear();
    flat_.resize(static_cast<size_t>(total_positions));

    // =========================================================================
    // Phase 2: re-scan packed genome and fill flat_ positions.
    // For each valid k-mer at position i:
    //   flat_[start + cursor] = i;  ++cursor;
    // After the scan, cursor == count for every entry.
    // =========================================================================
    std::cerr << "[kmer_index] Phase 2/2: filling position array...\n";
    {
        uint32_t n_count = 0;
        for (uint32_t i = 0; i < static_cast<uint32_t>(k_); ++i)
            if (is_n_fast(i)) ++n_count;

        uint64_t key = 0;
        for (uint32_t i = 0; i < static_cast<uint32_t>(k_); ++i)
            key = ((key << 2) | packed_base(i)) & kmer_mask;

        if (n_count == 0) {
            auto it = table_.find(key);
            if (it != table_.end())
                flat_[it->second.first + it->second.second++] = 0;
        }

        uint32_t next_report = REPORT_STEP;
        for (uint32_t i = 1; i <= last_start; ++i) {
            if (is_n_fast(i - 1))                             --n_count;
            if (is_n_fast(i + static_cast<uint32_t>(k_) - 1)) ++n_count;
            key = ((key << 2) | packed_base(i + static_cast<uint32_t>(k_) - 1))
                  & kmer_mask;

            if (n_count == 0) {
                auto it = table_.find(key);
                if (it != table_.end())
                    flat_[it->second.first + it->second.second++] = i;
            }

            if (i >= next_report) {
                std::cerr << "[kmer_index] phase 2: " << (i / 1'000'000u)
                          << " / " << (len / 1'000'000u) << " Mbp\n";
                next_report += REPORT_STEP;
            }
        }
    }

    std::cerr << "[kmer_index] k=" << k_
              << "  distinct_kmers=" << table_.size()
              << "  total_positions=" << flat_.size()
              << "  genome_size=" << len << "\n";
}

// ---------------------------------------------------------------------------
// save()
//
// Binary layout (all multi-byte values little-endian, mmap-ready):
//
//  [0..7]                    magic "TEB_IDX1"
//  [8]                       algo byte (0 = KMER)
//  [9]                       reserved  (0)
//  [10..17]                  genome_len  (uint64)
//  [18 .. 18+P-1]            packed genome  (P = (len+3)/4 bytes)
//  [pad to 8-byte boundary]  alignment for n_mask
//  [nmask_off .. +N-1]       n_mask  (N = ((len+63)/64)*8 bytes)
//  [+8]                      k  (uint64)
//  [+8]                      n_entries  (uint64)
//  [n_entries × 20 bytes]    bucket array: each entry is
//                              key(8) | n_pos(4) | ovfl_byte_offset(8)
//  [...]                     overflow array: uint32_t positions, all buckets
//                            concatenated in iteration order.
// ---------------------------------------------------------------------------
void KmerIndex::save(const std::string& idx_path) const {
    if (!packed_ptr_ || !nmask_ptr_)
        throw std::runtime_error("KmerIndex::save: index not built or loaded");

    const uint64_t gen_len      = genome_size_;
    const uint64_t packed_bytes = (gen_len + 3) / 4;
    const uint64_t nmask_bytes  = ((gen_len + 63) / 64) * 8;
    const uint64_t nmask_offset = ((18 + packed_bytes) + 7) & ~uint64_t{7};
    const uint64_t pad_bytes    = nmask_offset - (18 + packed_bytes);
    const uint64_t n_entries    = static_cast<uint64_t>(table_.size());

    std::ofstream out(idx_path, std::ios::binary);
    if (!out) throw std::runtime_error("KmerIndex::save: cannot open " + idx_path);

    // Magic
    out.write(MAGIC, 8);

    // Algo + reserved
    const uint8_t algo_byte = static_cast<uint8_t>(IndexAlgo::KMER);
    const uint8_t reserved  = 0;
    out.write(reinterpret_cast<const char*>(&algo_byte), 1);
    out.write(reinterpret_cast<const char*>(&reserved),  1);

    // genome_len
    out.write(reinterpret_cast<const char*>(&gen_len), 8);

    // packed genome
    out.write(reinterpret_cast<const char*>(packed_ptr_),
              static_cast<std::streamsize>(packed_bytes));

    // alignment padding
    static const uint8_t zeros[8] = {};
    if (pad_bytes > 0)
        out.write(reinterpret_cast<const char*>(zeros),
                  static_cast<std::streamsize>(pad_bytes));

    // n_mask
    out.write(reinterpret_cast<const char*>(nmask_ptr_),
              static_cast<std::streamsize>(nmask_bytes));

    // k
    const uint64_t k64 = static_cast<uint64_t>(k_);
    out.write(reinterpret_cast<const char*>(&k64), 8);

    // n_entries
    out.write(reinterpret_cast<const char*>(&n_entries), 8);

    // Bucket array (20 bytes each: key[8] n_pos[4] ovfl_off[8]).
    // The CSR layout stores table_[key] = {start_in_flat, count}.
    // ovfl_off = start * sizeof(uint32_t) (byte offset into flat_ array).
    for (const auto& [k, sc] : table_) {
        const uint32_t n_pos    = sc.second;
        const uint64_t ovfl_off = static_cast<uint64_t>(sc.first) * sizeof(uint32_t);
        out.write(reinterpret_cast<const char*>(&k),       8);
        out.write(reinterpret_cast<const char*>(&n_pos),   4);
        out.write(reinterpret_cast<const char*>(&ovfl_off),8);
    }

    // Overflow (position) array — single contiguous write.
    if (!flat_.empty())
        out.write(reinterpret_cast<const char*>(flat_.data()),
                  static_cast<std::streamsize>(flat_.size() * sizeof(uint32_t)));

    if (!out) throw std::runtime_error("KmerIndex::save: write error");
}

// ---------------------------------------------------------------------------
// load()
//
// mmap the index file.  SPEED mode prefaults all pages; MEMORY mode uses
// lazy demand-paging.  The packed genome and n_mask pointers are aimed
// directly at the mmap region (zero copies for genome data).  The hash
// table is rebuilt by iterating the serialised bucket array.
// ---------------------------------------------------------------------------
void KmerIndex::load(const std::string& idx_path, LoadMode mode) {
    // Release any previous mmap.
    if (mmap_ptr_ && mmap_ptr_ != MAP_FAILED) {
        munmap(mmap_ptr_, mmap_size_);
        mmap_ptr_ = nullptr;
    }
    if (mmap_fd_ >= 0) {
        close(mmap_fd_);
        mmap_fd_ = -1;
    }
    table_.clear();

    const int fd = open(idx_path.c_str(), O_RDONLY);
    if (fd < 0)
        throw std::runtime_error("KmerIndex::load: cannot open '" + idx_path
                                 + "': " + strerror(errno));

    struct stat st{};
    if (fstat(fd, &st) < 0) {
        close(fd);
        throw std::runtime_error("KmerIndex::load: fstat failed: "
                                 + std::string(strerror(errno)));
    }
    const size_t file_size = static_cast<size_t>(st.st_size);
    if (file_size < 18)  {
        close(fd);
        throw std::runtime_error("KmerIndex::load: file too small");
    }

    // Map the file.
    int mmap_flags = MAP_PRIVATE;
#ifdef MAP_POPULATE
    if (mode == LoadMode::SPEED) mmap_flags |= MAP_POPULATE;
#endif
    void* ptr = mmap(nullptr, file_size, PROT_READ, mmap_flags, fd, 0);
    if (ptr == MAP_FAILED) {
        close(fd);
        throw std::runtime_error("KmerIndex::load: mmap failed: "
                                 + std::string(strerror(errno)));
    }
#ifndef MAP_POPULATE
    // macOS / BSD: advise the kernel to prefault pages for SPEED mode.
    if (mode == LoadMode::SPEED)
        madvise(ptr, file_size, MADV_WILLNEED);
#endif

    mmap_ptr_  = ptr;
    mmap_size_ = file_size;
    mmap_fd_   = fd;

    const uint8_t* p   = static_cast<const uint8_t*>(ptr);
    size_t         off = 0;

    // Verify magic.
    if (std::memcmp(p + off, MAGIC, 8) != 0)
        throw std::runtime_error("KmerIndex::load: bad magic bytes (not a TEB_IDX1 file)");
    off += 8;

    // Verify algo byte.
    const uint8_t algo_byte = p[off++];
    if (algo_byte != static_cast<uint8_t>(IndexAlgo::KMER))
        throw std::runtime_error("KmerIndex::load: wrong algo byte — expected KMER");
    ++off; // skip reserved

    // Genome length.
    const uint64_t gen_len = load_le64(p + off); off += 8;
    genome_size_ = static_cast<uint32_t>(gen_len);

    const uint64_t packed_bytes = (gen_len + 3) / 4;
    const uint64_t nmask_offset = ((18 + packed_bytes) + 7) & ~uint64_t{7};
    const uint64_t nmask_bytes  = ((gen_len + 63) / 64) * 8;

    // Zero-copy genome pointers directly into the mmap region.
    packed_ptr_ = p + off;                                            // off = 18
    nmask_ptr_  = reinterpret_cast<const uint64_t*>(p + nmask_offset);
    off = nmask_offset + nmask_bytes;

    // k.
    k_ = static_cast<size_t>(load_le64(p + off)); off += 8;

    // n_entries.
    const uint64_t n_entries = load_le64(p + off); off += 8;

    // Bucket array starts at `off`; flat position array follows immediately.
    const uint8_t* buckets  = p + off;
    constexpr uint64_t BUCKET_STRIDE = 20; // key(8) + n_pos(4) + ovfl_off(8)
    const uint8_t* overflow = buckets + n_entries * BUCKET_STRIDE;

    // Reconstruct the flat_ array (single memcpy from mmap region).
    const size_t ovfl_bytes = file_size - static_cast<size_t>(overflow - p);
    flat_.resize(ovfl_bytes / sizeof(uint32_t));
    if (!flat_.empty())
        std::memcpy(flat_.data(), overflow, ovfl_bytes);

    // Rebuild the hash table: table_[key] = {start_in_flat, count}.
    table_.reserve(n_entries);
    for (uint64_t i = 0; i < n_entries; ++i) {
        const uint8_t* be  = buckets + i * BUCKET_STRIDE;
        const uint64_t key      = load_le64(be);
        const uint32_t n_pos    = load_le32(be + 8);
        const uint64_t ovfl_off = load_le64(be + 12);
        const uint32_t start    = static_cast<uint32_t>(ovfl_off / sizeof(uint32_t));
        table_.emplace(key, std::make_pair(start, n_pos));
    }
}

// ---------------------------------------------------------------------------
// query()
// ---------------------------------------------------------------------------
std::vector<uint32_t> KmerIndex::query(const std::string& pattern) const {
    if (pattern.size() != k_) return {};

    const uint64_t kmer_mask = (k_ < 32)
        ? ((uint64_t{1} << (2 * k_)) - 1)
        : ~uint64_t{0};

    uint64_t key = 0;
    for (char c : pattern) {
        const uint8_t b = BASE_ENC.v[static_cast<uint8_t>(std::toupper(c))];
        if (b == 0xFFu) return {}; // N or unknown base
        key = ((key << 2) | b) & kmer_mask;
    }

    const auto it = table_.find(key);
    if (it == table_.end()) return {};
    const auto [start, cnt] = it->second;
    if (cnt == 0) return {};
    return std::vector<uint32_t>(flat_.begin() + start,
                                 flat_.begin() + start + cnt);
}

// ---------------------------------------------------------------------------
// verify_match() — delegates to count_mismatches() from sequence.hpp
// ---------------------------------------------------------------------------
int KmerIndex::verify_match(uint32_t genome_pos,
                             const uint8_t* read_packed,
                             uint32_t read_len) const {
    return count_mismatches(packed_ptr_, nmask_ptr_,
                            genome_pos, read_packed, read_len);
}


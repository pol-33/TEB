#include "seed_index.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

constexpr char kMagic[8] = {'K', '1', '4', 'I', 'D', 'X', '\0', '\0'};
constexpr uint32_t kVersion = 1;
constexpr uint64_t kSeedMask = kSeedBucketCount - 1ULL;
constexpr uint32_t kPositionWritePartitions = 8;

struct WritableMapping {
    void* base = MAP_FAILED;
    std::size_t mapped_size = 0;
    uint8_t* data = nullptr;
};

void append_bytes(std::vector<uint8_t>& buffer, const void* data, std::size_t size) {
    const auto* bytes = static_cast<const uint8_t*>(data);
    buffer.insert(buffer.end(), bytes, bytes + size);
}

template <typename T>
void append_value(std::vector<uint8_t>& buffer, const T& value) {
    append_bytes(buffer, &value, sizeof(T));
}

template <typename T>
T read_value(const uint8_t*& ptr, const uint8_t* end) {
    if (static_cast<std::size_t>(end - ptr) < sizeof(T)) {
        throw std::runtime_error("corrupt seed index");
    }
    T value{};
    std::memcpy(&value, ptr, sizeof(T));
    ptr += sizeof(T);
    return value;
}

std::string normalize_name(std::string_view raw_name) {
    const std::size_t cut = raw_name.find_first_of(" \t");
    const std::string_view trimmed = (cut == std::string_view::npos) ? raw_name : raw_name.substr(0, cut);
    return std::string(trimmed);
}

uint8_t two_bit_base(char base) {
    return packed_base_from_code(encode_base(base));
}

void set_packed_code_raw(uint8_t* packed, uint64_t index, uint8_t code) {
    const std::size_t byte_index = static_cast<std::size_t>(index / 4ULL);
    const uint8_t shift = static_cast<uint8_t>((index % 4ULL) * 2ULL);
    const uint8_t value = packed_base_from_code(code);
    packed[byte_index] &= static_cast<uint8_t>(~(0x3u << shift));
    packed[byte_index] |= static_cast<uint8_t>(value << shift);
}

void write_all_at(int fd, const void* data, std::size_t bytes, uint64_t offset) {
    const auto* ptr = static_cast<const uint8_t*>(data);
    std::size_t written = 0;
    while (written < bytes) {
        const ssize_t rc = ::pwrite(fd, ptr + written, bytes - written, static_cast<off_t>(offset + written));
        if (rc < 0) {
            throw std::system_error(errno, std::generic_category(), "pwrite failed while building seed index");
        }
        written += static_cast<std::size_t>(rc);
    }
}

WritableMapping map_writable_region(int fd, uint64_t offset, std::size_t size) {
    if (size == 0) {
        return WritableMapping{};
    }

    const long page_size_long = ::sysconf(_SC_PAGESIZE);
    if (page_size_long <= 0) {
        throw std::runtime_error("failed to determine system page size");
    }

    const uint64_t page_size = static_cast<uint64_t>(page_size_long);
    const uint64_t aligned_offset = (offset / page_size) * page_size;
    const std::size_t delta = static_cast<std::size_t>(offset - aligned_offset);
    const std::size_t mapped_size = delta + size;

    void* base = ::mmap(nullptr,
                        mapped_size,
                        PROT_READ | PROT_WRITE,
                        MAP_SHARED,
                        fd,
                        static_cast<off_t>(aligned_offset));
    if (base == MAP_FAILED) {
        throw std::system_error(errno, std::generic_category(), "mmap failed while building seed index");
    }

    WritableMapping mapping;
    mapping.base = base;
    mapping.mapped_size = mapped_size;
    mapping.data = static_cast<uint8_t*>(base) + delta;
    return mapping;
}

void unmap_writable_region(WritableMapping& mapping) noexcept {
    if (mapping.base != MAP_FAILED) {
        ::munmap(mapping.base, mapping.mapped_size);
        mapping.base = MAP_FAILED;
        mapping.mapped_size = 0;
        mapping.data = nullptr;
    }
}

struct PassOneResult {
    std::vector<ChromInfo> chromosomes;
    uint64_t genome_length = 0;
};

PassOneResult scan_reference_and_count(const std::string& fasta_path, std::vector<uint32_t>& offsets) {
    std::ifstream in(fasta_path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA file: " + fasta_path);
    }

    PassOneResult result;
    std::string line;
    ChromInfo current;
    bool have_current = false;
    uint64_t rolling = 0;
    uint32_t valid = 0;
    uint64_t next_report = 250000000ULL;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            if (have_current) {
                result.chromosomes.push_back(current);
            }
            current = ChromInfo{};
            current.name = normalize_name(line.substr(1));
            current.offset = result.genome_length;
            current.length = 0;
            have_current = true;
            rolling = 0;
            valid = 0;
            continue;
        }
        if (!have_current) {
            throw std::runtime_error("FASTA sequence found before any header in: " + fasta_path);
        }

        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
                continue;
            }
            const uint8_t base = two_bit_base(ch);
            rolling = ((rolling << 2U) | static_cast<uint64_t>(base)) & kSeedMask;
            ++valid;
            ++current.length;
            ++result.genome_length;
            if (valid >= kSeedLength) {
                ++offsets[static_cast<std::size_t>(rolling)];
            }
            if (result.genome_length >= next_report) {
                std::cerr << "[indexer] counted " << result.genome_length << " bases\n";
                next_report += 250000000ULL;
            }
        }
    }

    if (have_current) {
        result.chromosomes.push_back(current);
    }
    if (result.chromosomes.empty()) {
        throw std::runtime_error("no chromosomes found in FASTA file: " + fasta_path);
    }

    return result;
}

void populate_packed_genome(const std::string& fasta_path, uint8_t* packed_genome) {
    std::ifstream in(fasta_path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA file while writing packed genome: " + fasta_path);
    }

    std::string line;
    bool have_current = false;
    uint64_t genome_pos = 0;
    uint64_t next_report = 250000000ULL;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            have_current = true;
            continue;
        }
        if (!have_current) {
            throw std::runtime_error("FASTA sequence found before any header while writing packed genome");
        }

        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
                continue;
            }
            const uint8_t code = encode_base(ch);
            set_packed_code_raw(packed_genome, genome_pos, code);

            ++genome_pos;
            if (genome_pos >= next_report) {
                std::cerr << "[indexer] wrote " << genome_pos << " packed reference bases\n";
                std::cerr.flush();
                next_report += 250000000ULL;
            }
        }
    }
}

void populate_positions_partition(const std::string& fasta_path,
                                  uint32_t partition_begin_code,
                                  uint32_t partition_end_code,
                                  uint32_t* positions,
                                  std::vector<uint32_t>& cursors,
                                  uint32_t partition_index) {
    std::ifstream in(fasta_path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA file while writing position partition: " + fasta_path);
    }

    std::string line;
    bool have_current = false;
    uint64_t rolling = 0;
    uint32_t valid = 0;
    uint64_t genome_pos = 0;
    uint64_t next_report = 500000000ULL;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            have_current = true;
            rolling = 0;
            valid = 0;
            continue;
        }
        if (!have_current) {
            throw std::runtime_error("FASTA sequence found before any header while writing position partition");
        }

        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
                continue;
            }

            const uint8_t code = encode_base(ch);
            rolling = ((rolling << 2U) | static_cast<uint64_t>(packed_base_from_code(code))) & kSeedMask;
            ++valid;
            if (valid >= kSeedLength &&
                rolling >= static_cast<uint64_t>(partition_begin_code) &&
                rolling < static_cast<uint64_t>(partition_end_code)) {
                const uint32_t local_code = static_cast<uint32_t>(rolling) - partition_begin_code;
                const uint32_t slot = cursors[static_cast<std::size_t>(local_code)]++;
                positions[slot] = static_cast<uint32_t>(genome_pos - static_cast<uint64_t>(kSeedLength) + 1ULL);
            }

            ++genome_pos;
            if (genome_pos >= next_report) {
                std::cerr << "[indexer] partition " << (partition_index + 1U) << "/" << kPositionWritePartitions
                          << " scanned " << genome_pos << " bases\n";
                std::cerr.flush();
                next_report += 500000000ULL;
            }
        }
    }
}

}  // namespace

uint32_t encode_seed(std::string_view seed) {
    if (seed.size() != kSeedLength) {
        throw std::runtime_error("seed length mismatch");
    }
    uint32_t code = 0;
    for (char base : seed) {
        code = static_cast<uint32_t>((code << 2U) | two_bit_base(base));
    }
    return code;
}

IndexBuildStats build_seed_index(const std::string& fasta_path, const std::string& index_path) {
    std::cerr << "[indexer] counting " << kSeedLength << "-mer seeds across the reference\n";
    std::vector<uint32_t> offsets(static_cast<std::size_t>(kSeedBucketCount) + 1U, 0U);
    PassOneResult pass_one = scan_reference_and_count(fasta_path, offsets);

    uint64_t indexed_positions = 0;
    for (uint64_t bucket = 0; bucket < kSeedBucketCount; ++bucket) {
        const uint32_t count = offsets[static_cast<std::size_t>(bucket)];
        offsets[static_cast<std::size_t>(bucket)] = static_cast<uint32_t>(indexed_positions);
        indexed_positions += count;
    }
    offsets[static_cast<std::size_t>(kSeedBucketCount)] = static_cast<uint32_t>(indexed_positions);

    std::cerr << "[indexer] counted " << indexed_positions << " indexed seed positions\n";

    std::vector<uint8_t> header;
    append_bytes(header, kMagic, sizeof(kMagic));
    append_value(header, kVersion);
    append_value(header, kSeedLength);
    append_value(header, pass_one.genome_length);
    append_value(header, indexed_positions);
    const uint32_t num_chromosomes = static_cast<uint32_t>(pass_one.chromosomes.size());
    append_value(header, num_chromosomes);
    for (const ChromInfo& chrom : pass_one.chromosomes) {
        if (chrom.name.size() > 0xFFFFu) {
            throw std::runtime_error("chromosome name too long for serialized index");
        }
        const uint16_t name_len = static_cast<uint16_t>(chrom.name.size());
        append_value(header, name_len);
        append_bytes(header, chrom.name.data(), chrom.name.size());
        append_value(header, chrom.offset);
        append_value(header, chrom.length);
    }

    const uint64_t packed_genome_bytes = packed_byte_count(pass_one.genome_length);
    const uint64_t offsets_bytes = (kSeedBucketCount + 1ULL) * sizeof(uint32_t);
    const uint64_t positions_bytes = indexed_positions * sizeof(uint32_t);
    const uint64_t packed_genome_offset = header.size();
    const uint64_t offsets_offset = packed_genome_offset + packed_genome_bytes;
    const uint64_t positions_offset = offsets_offset + offsets_bytes;
    const uint64_t file_size = positions_offset + positions_bytes;

    const int fd = ::open(index_path.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fd < 0) {
        throw std::system_error(errno, std::generic_category(), "failed to open output seed index");
    }

    WritableMapping packed_genome_map;
    WritableMapping positions_map;
    try {
#ifdef F_NOCACHE
        (void)::fcntl(fd, F_NOCACHE, 1);
#endif
        if (::ftruncate(fd, static_cast<off_t>(file_size)) != 0) {
            throw std::system_error(errno, std::generic_category(), "ftruncate failed while sizing seed index");
        }

        write_all_at(fd, header.data(), header.size(), 0);
        write_all_at(fd, offsets.data(), static_cast<std::size_t>(offsets_bytes), offsets_offset);

        packed_genome_map = map_writable_region(fd, packed_genome_offset, static_cast<std::size_t>(packed_genome_bytes));
        std::cerr << "[indexer] populating packed genome\n";
        populate_packed_genome(fasta_path, packed_genome_map.data);
        unmap_writable_region(packed_genome_map);

        std::cerr << "[indexer] populating on-disk seed positions across "
                  << kPositionWritePartitions << " partition(s)\n";
        const uint64_t codes_per_partition = kSeedBucketCount / static_cast<uint64_t>(kPositionWritePartitions);
        for (uint32_t partition = 0; partition < kPositionWritePartitions; ++partition) {
            const uint32_t begin_code = static_cast<uint32_t>(static_cast<uint64_t>(partition) * codes_per_partition);
            const uint32_t end_code = (partition + 1U == kPositionWritePartitions)
                ? static_cast<uint32_t>(kSeedBucketCount)
                : static_cast<uint32_t>(static_cast<uint64_t>(partition + 1U) * codes_per_partition);
            const uint32_t begin_slot = offsets[begin_code];
            const uint32_t end_slot = offsets[end_code];
            const uint64_t partition_positions = static_cast<uint64_t>(end_slot) - static_cast<uint64_t>(begin_slot);

            std::cerr << "[indexer] writing partition " << (partition + 1U) << "/" << kPositionWritePartitions
                      << " with " << partition_positions << " seed hits\n";
            std::cerr.flush();

            if (partition_positions == 0) {
                continue;
            }

            positions_map = map_writable_region(fd,
                                                positions_offset + static_cast<uint64_t>(begin_slot) * sizeof(uint32_t),
                                                static_cast<std::size_t>(partition_positions * sizeof(uint32_t)));

            std::vector<uint32_t> partition_cursors(static_cast<std::size_t>(end_code - begin_code), 0U);
            for (uint32_t code = begin_code; code < end_code; ++code) {
                partition_cursors[static_cast<std::size_t>(code - begin_code)] = offsets[code] - begin_slot;
            }

            populate_positions_partition(fasta_path,
                                         begin_code,
                                         end_code,
                                         reinterpret_cast<uint32_t*>(positions_map.data),
                                         partition_cursors,
                                         partition);
            unmap_writable_region(positions_map);
        }
        ::close(fd);
    } catch (...) {
        unmap_writable_region(packed_genome_map);
        unmap_writable_region(positions_map);
        ::close(fd);
        throw;
    }

    IndexBuildStats stats;
    stats.genome_length = pass_one.genome_length;
    stats.indexed_positions = indexed_positions;
    stats.num_chromosomes = num_chromosomes;
    stats.file_size = file_size;
    return stats;
}

SeedIndexView::SeedIndexView()
    : fd_(-1),
      mapping_(nullptr),
      file_size_(0),
      seed_length_(0),
      genome_length_(0),
      indexed_positions_(0),
      packed_genome_(nullptr),
      offsets_(nullptr),
      positions_(nullptr) {}

SeedIndexView::SeedIndexView(const std::string& path) : SeedIndexView() {
    open(path);
}

SeedIndexView::~SeedIndexView() {
    close();
}

SeedIndexView::SeedIndexView(SeedIndexView&& other) noexcept : SeedIndexView() {
    move_from(std::move(other));
}

SeedIndexView& SeedIndexView::operator=(SeedIndexView&& other) noexcept {
    if (this != &other) {
        close();
        move_from(std::move(other));
    }
    return *this;
}

void SeedIndexView::move_from(SeedIndexView&& other) noexcept {
    fd_ = other.fd_;
    mapping_ = other.mapping_;
    file_size_ = other.file_size_;
    seed_length_ = other.seed_length_;
    genome_length_ = other.genome_length_;
    indexed_positions_ = other.indexed_positions_;
    packed_genome_ = other.packed_genome_;
    offsets_ = other.offsets_;
    positions_ = other.positions_;
    chromosomes_ = std::move(other.chromosomes_);

    other.fd_ = -1;
    other.mapping_ = nullptr;
    other.file_size_ = 0;
    other.seed_length_ = 0;
    other.genome_length_ = 0;
    other.indexed_positions_ = 0;
    other.packed_genome_ = nullptr;
    other.offsets_ = nullptr;
    other.positions_ = nullptr;
    other.chromosomes_.clear();
}

void SeedIndexView::open(const std::string& path) {
    close();

    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("failed to open seed index: " + path);
    }

    struct stat st {};
    if (::fstat(fd_, &st) != 0) {
        const int err = errno;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "fstat failed for seed index");
    }

    file_size_ = static_cast<std::size_t>(st.st_size);
    mapping_ = static_cast<const uint8_t*>(::mmap(nullptr, file_size_, PROT_READ, MAP_SHARED, fd_, 0));
    if (mapping_ == MAP_FAILED) {
        const int err = errno;
        mapping_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "mmap failed for seed index");
    }

#ifdef MADV_RANDOM
    (void)::madvise(const_cast<uint8_t*>(mapping_), file_size_, MADV_RANDOM);
#endif
    parse();
}

void SeedIndexView::close() {
    chromosomes_.clear();
    if (mapping_ != nullptr) {
        ::munmap(const_cast<uint8_t*>(mapping_), file_size_);
        mapping_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
    file_size_ = 0;
    seed_length_ = 0;
    genome_length_ = 0;
    indexed_positions_ = 0;
    packed_genome_ = nullptr;
    offsets_ = nullptr;
    positions_ = nullptr;
}

bool SeedIndexView::is_open() const {
    return mapping_ != nullptr;
}

void SeedIndexView::parse() {
    const uint8_t* ptr = mapping_;
    const uint8_t* end = mapping_ + file_size_;

    if (static_cast<std::size_t>(end - ptr) < sizeof(kMagic) || std::memcmp(ptr, kMagic, sizeof(kMagic)) != 0) {
        throw std::runtime_error("invalid seed index magic header");
    }
    ptr += sizeof(kMagic);

    const uint32_t version = read_value<uint32_t>(ptr, end);
    if (version != kVersion) {
        throw std::runtime_error("unsupported seed index version");
    }

    seed_length_ = read_value<uint32_t>(ptr, end);
    genome_length_ = read_value<uint64_t>(ptr, end);
    indexed_positions_ = read_value<uint64_t>(ptr, end);
    const uint32_t num_chromosomes = read_value<uint32_t>(ptr, end);

    chromosomes_.clear();
    chromosomes_.reserve(num_chromosomes);
    for (uint32_t i = 0; i < num_chromosomes; ++i) {
        const uint16_t name_len = read_value<uint16_t>(ptr, end);
        if (static_cast<std::size_t>(end - ptr) < name_len) {
            throw std::runtime_error("corrupt chromosome table in seed index");
        }
        const char* name_ptr = reinterpret_cast<const char*>(ptr);
        ptr += name_len;
        ChromView chrom;
        chrom.name = std::string_view(name_ptr, name_len);
        chrom.offset = read_value<uint64_t>(ptr, end);
        chrom.length = read_value<uint64_t>(ptr, end);
        chromosomes_.push_back(chrom);
    }

    const std::size_t packed_genome_bytes = packed_byte_count(genome_length_);
    const std::size_t offsets_bytes = static_cast<std::size_t>((kSeedBucketCount + 1ULL) * sizeof(uint32_t));
    const std::size_t positions_bytes = static_cast<std::size_t>(indexed_positions_ * sizeof(uint32_t));
    if (static_cast<std::size_t>(end - ptr) < packed_genome_bytes + offsets_bytes + positions_bytes) {
        throw std::runtime_error("seed index truncated");
    }

    packed_genome_ = ptr;
    ptr += packed_genome_bytes;
    offsets_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += offsets_bytes;
    positions_ = reinterpret_cast<const uint32_t*>(ptr);
}

uint32_t SeedIndexView::seed_length() const {
    return seed_length_;
}

uint64_t SeedIndexView::genome_length() const {
    return genome_length_;
}

uint64_t SeedIndexView::indexed_positions() const {
    return indexed_positions_;
}

const std::vector<SeedIndexView::ChromView>& SeedIndexView::chromosomes() const {
    return chromosomes_;
}

uint32_t SeedIndexView::bucket_begin(uint32_t code) const {
    return offsets_[code];
}

uint32_t SeedIndexView::bucket_end(uint32_t code) const {
    return offsets_[static_cast<std::size_t>(code) + 1U];
}

uint32_t SeedIndexView::position_at(uint64_t idx) const {
    return positions_[idx];
}

uint64_t SeedIndexView::chromosome_index(uint64_t genome_pos) const {
    auto it = std::upper_bound(chromosomes_.begin(), chromosomes_.end(), genome_pos, [](uint64_t pos, const ChromView& chrom) {
        return pos < chrom.offset;
    });
    if (it == chromosomes_.begin()) {
        throw std::runtime_error("genome position before first chromosome");
    }
    --it;
    if (genome_pos >= it->offset + it->length) {
        throw std::runtime_error("genome position outside chromosome bounds");
    }
    return static_cast<uint64_t>(std::distance(chromosomes_.begin(), it));
}

bool SeedIndexView::same_chromosome_window(uint64_t genome_pos, uint64_t span) const {
    if (span == 0 || genome_pos >= genome_length_) {
        return false;
    }
    const uint64_t idx = chromosome_index(genome_pos);
    const ChromView& chrom = chromosomes_[static_cast<std::size_t>(idx)];
    return genome_pos + span <= chrom.offset + chrom.length;
}

std::string_view SeedIndexView::chromosome_name(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].name;
}

uint64_t SeedIndexView::chromosome_offset(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].offset;
}

uint64_t SeedIndexView::chromosome_length(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].length;
}

char SeedIndexView::reference_base(uint64_t genome_pos) const {
    return decode_base(get_packed_code(packed_genome_, genome_pos));
}

void SeedIndexView::extract_reference(uint64_t genome_pos, uint64_t length, std::string& out) const {
    out.resize(static_cast<std::size_t>(length));
    for (uint64_t i = 0; i < length; ++i) {
        out[static_cast<std::size_t>(i)] = reference_base(genome_pos + i);
    }
}

}  // namespace mapper_memory

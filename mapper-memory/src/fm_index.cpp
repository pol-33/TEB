#include "fm_index.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <system_error>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "nucleotide.hpp"
#include "simd_dispatch.hpp"

namespace mapper_memory {

namespace {

// Updated magic with version for new format
constexpr char kMagic[8] = {'F', 'M', 'C', 'H', 'R', '0', '0', '2'};  // Version 002

// Feature flags
constexpr uint8_t kFlagHasGenome = 0x01;

struct ChromosomeHeader {
    uint16_t name_len = 0;
    uint8_t flags = 0;          // Per-chromosome flags
    uint8_t reserved = 0;
    uint32_t reserved2 = 0;
    uint64_t length = 0;
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    uint64_t packed_bwt_offset = 0;
    uint64_t packed_bwt_bytes = 0;
    uint64_t occ_offset = 0;
    uint64_t occ_entries = 0;
    uint64_t sa_offset = 0;
    uint64_t sa_entries = 0;
    uint64_t genome_offset = 0;
    uint64_t genome_bytes = 0;
    uint64_t c_array[kAlphabetSize]{};
};

template <typename T>
void write_value(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) {
        throw std::runtime_error("failed while writing FM-index file");
    }
}

template <typename T>
T read_value(const uint8_t*& ptr, const uint8_t* end) {
    if (static_cast<std::size_t>(end - ptr) < sizeof(T)) {
        throw std::runtime_error("corrupt FM-index file");
    }
    T value{};
    std::memcpy(&value, ptr, sizeof(T));
    ptr += sizeof(T);
    return value;
}

std::size_t occ_checkpoint_count(uint64_t n, uint32_t occ_sample) {
    return static_cast<std::size_t>(n / occ_sample) + 1U;
}

std::size_t sa_checkpoint_count(uint64_t n, uint32_t sa_sample) {
    return static_cast<std::size_t>(n / sa_sample) + 1U;
}

uint32_t dna_occ_index(uint8_t rank) {
    return static_cast<uint32_t>(rank - kARank);
}

std::vector<uint8_t> pack_sequence(const std::string& sequence) {
    std::vector<uint8_t> packed;
    packed.reserve(packed_byte_count(sequence.size()));
    for (uint64_t i = 0; i < sequence.size(); ++i) {
        append_packed_rank(packed, i, rank_from_base(sequence[static_cast<std::size_t>(i)]));
    }
    return packed;
}

uint64_t align_up(uint64_t value, uint64_t alignment) {
    return (value + alignment - 1U) & ~(alignment - 1U);
}

}  // namespace

OwnedChromosomeIndex build_chromosome_index(const std::string& name,
                                            const std::string& sequence,
                                            const BWTData& bwt,
                                            const std::vector<uint32_t>& suffix_array,
                                            const IndexConfig& config) {
    OwnedChromosomeIndex index;
    index.name = name;
    index.length = sequence.size();
    index.text_length = bwt.text_length;
    index.primary_index = bwt.primary_index;
    index.occ_sample = config.occ_sample;
    index.sa_sample = config.sa_sample;
    index.packed_bwt = bwt.packed_bwt;
    
    // Only pack genome if configured
    if (config.store_genome) {
        index.packed_genome = pack_sequence(sequence);
    }

    // Build C array
    index.c_array[0] = 0;
    for (std::size_t rank = 1; rank < index.c_array.size(); ++rank) {
        index.c_array[rank] = index.c_array[rank - 1U] + bwt.counts[rank - 1U];
    }

    // Build Occ checkpoints
    const std::size_t occ_blocks = occ_checkpoint_count(index.text_length, config.occ_sample);
    index.sampled_occ.assign(occ_blocks * 4U, 0U);
    std::array<uint32_t, 4> running{};
    for (uint64_t row = 0; row < index.text_length; ++row) {
        const uint8_t rank = (row == index.primary_index) ? kSentinelRank : get_packed_rank(index.packed_bwt.data(), row);
        if (rank >= kARank) {
            ++running[dna_occ_index(rank)];
        }
        if ((row + 1ULL) % config.occ_sample == 0ULL) {
            const std::size_t block = static_cast<std::size_t>((row + 1ULL) / config.occ_sample);
            for (std::size_t i = 0; i < running.size(); ++i) {
                index.sampled_occ[block * 4U + i] = running[i];
            }
        }
    }

    // Build SA samples
    const std::size_t sa_entries = sa_checkpoint_count(index.text_length, config.sa_sample);
    index.sampled_sa.assign(sa_entries, 0U);
    for (uint64_t row = 0; row < index.text_length; row += config.sa_sample) {
        index.sampled_sa[static_cast<std::size_t>(row / config.sa_sample)] = suffix_array[static_cast<std::size_t>(row)];
    }
    
    return index;
}

OwnedChromosomeIndex build_chromosome_index(const std::string& name,
                                            const std::string& sequence,
                                            const BWTData& bwt,
                                            const std::vector<uint32_t>& suffix_array,
                                            uint32_t occ_sample,
                                            uint32_t sa_sample) {
    IndexConfig config;
    config.occ_sample = occ_sample;
    config.sa_sample = sa_sample;
    config.store_genome = true;
    return build_chromosome_index(name, sequence, bwt, suffix_array, config);
}

std::size_t FMIndexView::write(const std::string& path,
                               const std::vector<OwnedChromosomeIndex>& chromosomes,
                               const IndexConfig& config) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("failed to open FM-index output file: " + path);
    }

    // Build global flags
    uint8_t global_flags = 0;
    if (config.store_genome) global_flags |= kFlagHasGenome;

    const uint32_t num_chroms = static_cast<uint32_t>(chromosomes.size());
    const uint64_t prelude_bytes = sizeof(kMagic) + sizeof(uint32_t) * 3U + sizeof(uint8_t) * 4U;
    uint64_t header_bytes = prelude_bytes;
    for (const OwnedChromosomeIndex& chrom : chromosomes) {
        header_bytes += sizeof(ChromosomeHeader) + chrom.name.size();
    }
    const uint64_t data_start = align_up(header_bytes, 8U);

    std::vector<ChromosomeHeader> headers(chromosomes.size());
    uint64_t offset = data_start;
    for (std::size_t i = 0; i < chromosomes.size(); ++i) {
        const OwnedChromosomeIndex& chrom = chromosomes[i];
        ChromosomeHeader& header = headers[i];
        if (chrom.name.size() > 0xFFFFu) {
            throw std::runtime_error("chromosome name too long for serialized FM-index");
        }

        // Set per-chromosome flags
        header.flags = 0;
        if (!chrom.packed_genome.empty()) header.flags |= kFlagHasGenome;

        header.name_len = static_cast<uint16_t>(chrom.name.size());
        header.length = chrom.length;
        header.text_length = chrom.text_length;
        header.primary_index = chrom.primary_index;
        std::copy(chrom.c_array.begin(), chrom.c_array.end(), header.c_array);

        // Packed BWT
        header.packed_bwt_offset = offset;
        header.packed_bwt_bytes = chrom.packed_bwt.size();
        offset += chrom.packed_bwt.size();
        offset = align_up(offset, 8U);

        // Occ checkpoints
        header.occ_offset = offset;
        header.occ_entries = chrom.sampled_occ.size();
        offset += chrom.sampled_occ.size() * sizeof(uint32_t);
        offset = align_up(offset, 8U);

        // 32-bit SA
        header.sa_offset = offset;
        header.sa_entries = chrom.sampled_sa.size();
        offset += chrom.sampled_sa.size() * sizeof(uint32_t);
        offset = align_up(offset, 8U);

        // Packed genome (optional)
        header.genome_offset = offset;
        header.genome_bytes = chrom.packed_genome.size();
        offset += chrom.packed_genome.size();
        offset = align_up(offset, 8U);
    }

    // Write header
    out.write(kMagic, sizeof(kMagic));
    write_value(out, config.occ_sample);
    write_value(out, config.sa_sample);
    write_value(out, num_chroms);
    write_value(out, global_flags);
    uint8_t reserved[3] = {0, 0, 0};
    out.write(reinterpret_cast<const char*>(reserved), 3);
    
    for (std::size_t i = 0; i < chromosomes.size(); ++i) {
        out.write(reinterpret_cast<const char*>(&headers[i]), sizeof(ChromosomeHeader));
        out.write(chromosomes[i].name.data(), static_cast<std::streamsize>(chromosomes[i].name.size()));
    }

    // Padding to data_start
    const uint64_t current = static_cast<uint64_t>(out.tellp());
    if (current < data_start) {
        std::vector<char> padding(static_cast<std::size_t>(data_start - current), 0);
        out.write(padding.data(), static_cast<std::streamsize>(padding.size()));
    }

    auto write_aligned_blob = [&](const auto* data, std::size_t bytes) {
        if (bytes > 0U) {
            out.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(bytes));
        }
        const uint64_t pos = static_cast<uint64_t>(out.tellp());
        const uint64_t aligned = align_up(pos, 8U);
        if (aligned > pos) {
            std::vector<char> padding(static_cast<std::size_t>(aligned - pos), 0);
            out.write(padding.data(), static_cast<std::streamsize>(padding.size()));
        }
    };

    for (const OwnedChromosomeIndex& chrom : chromosomes) {
        write_aligned_blob(chrom.packed_bwt.data(), chrom.packed_bwt.size());
        write_aligned_blob(chrom.sampled_occ.data(), chrom.sampled_occ.size() * sizeof(uint32_t));
        write_aligned_blob(chrom.sampled_sa.data(), chrom.sampled_sa.size() * sizeof(uint32_t));
        write_aligned_blob(chrom.packed_genome.data(), chrom.packed_genome.size());
    }

    if (!out) {
        throw std::runtime_error("failed while writing FM-index payload");
    }
    return static_cast<std::size_t>(out.tellp());
}

std::size_t FMIndexView::write(const std::string& path,
                               const std::vector<OwnedChromosomeIndex>& chromosomes,
                               uint32_t occ_sample,
                               uint32_t sa_sample) {
    IndexConfig config;
    config.occ_sample = occ_sample;
    config.sa_sample = sa_sample;
    config.store_genome = true;
    return write(path, chromosomes, config);
}

FMIndexView::FMIndexView()
    : fd_(-1),
      mapping_(nullptr),
      file_size_(0),
      occ_sample_(0),
      sa_sample_(0),
      flags_(0) {}

FMIndexView::FMIndexView(const std::string& path) : FMIndexView() {
    open(path);
}

FMIndexView::~FMIndexView() {
    close();
}

FMIndexView::FMIndexView(FMIndexView&& other) noexcept : FMIndexView() {
    move_from(std::move(other));
}

FMIndexView& FMIndexView::operator=(FMIndexView&& other) noexcept {
    if (this != &other) {
        close();
        move_from(std::move(other));
    }
    return *this;
}

void FMIndexView::move_from(FMIndexView&& other) noexcept {
    fd_ = other.fd_;
    mapping_ = other.mapping_;
    file_size_ = other.file_size_;
    occ_sample_ = other.occ_sample_;
    sa_sample_ = other.sa_sample_;
    flags_ = other.flags_;
    chromosomes_ = std::move(other.chromosomes_);

    other.fd_ = -1;
    other.mapping_ = nullptr;
    other.file_size_ = 0;
    other.occ_sample_ = 0;
    other.sa_sample_ = 0;
    other.flags_ = 0;
    other.chromosomes_.clear();
}

void FMIndexView::open(const std::string& path) {
    close();

    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("failed to open FM-index: " + path);
    }

    struct stat st {};
    if (::fstat(fd_, &st) != 0) {
        const int err = errno;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "fstat failed for FM-index");
    }

    file_size_ = static_cast<std::size_t>(st.st_size);
    mapping_ = static_cast<const uint8_t*>(::mmap(nullptr, file_size_, PROT_READ, MAP_SHARED, fd_, 0));
    if (mapping_ == MAP_FAILED) {
        const int err = errno;
        mapping_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "mmap failed for FM-index");
    }

#ifdef MADV_RANDOM
    (void)::madvise(const_cast<uint8_t*>(mapping_), file_size_, MADV_RANDOM);
#endif
    parse();
}

void FMIndexView::close() {
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
    occ_sample_ = 0;
    sa_sample_ = 0;
    flags_ = 0;
}

void FMIndexView::parse() {
    const uint8_t* ptr = mapping_;
    const uint8_t* end = mapping_ + file_size_;

    // Check magic - support both old and new format
    if (static_cast<std::size_t>(end - ptr) < sizeof(kMagic)) {
        throw std::runtime_error("invalid FM-index magic header");
    }
    
    // Check for version 001 (old format) or 002 (new format)
    bool is_v1 = (std::memcmp(ptr, "FMCHR001", 8) == 0);
    bool is_v2 = (std::memcmp(ptr, kMagic, sizeof(kMagic)) == 0);
    if (!is_v1 && !is_v2) {
        throw std::runtime_error("invalid FM-index magic header");
    }
    ptr += sizeof(kMagic);

    occ_sample_ = read_value<uint32_t>(ptr, end);
    sa_sample_ = read_value<uint32_t>(ptr, end);
    const uint32_t num_chroms = read_value<uint32_t>(ptr, end);
    
    if (is_v2) {
        flags_ = read_value<uint8_t>(ptr, end);
        ptr += 3;  // Skip reserved bytes
    } else {
        flags_ = kFlagHasGenome;  // V1 always has genome
        (void)read_value<uint32_t>(ptr, end);  // Skip old reserved field
    }

    chromosomes_.clear();
    chromosomes_.reserve(num_chroms);
    for (uint32_t i = 0; i < num_chroms; ++i) {
        if (static_cast<std::size_t>(end - ptr) < sizeof(ChromosomeHeader)) {
            throw std::runtime_error("corrupt FM-index chromosome header");
        }
        ChromosomeHeader header{};
        std::memcpy(&header, ptr, sizeof(ChromosomeHeader));
        ptr += sizeof(ChromosomeHeader);
        if (static_cast<std::size_t>(end - ptr) < header.name_len) {
            throw std::runtime_error("corrupt FM-index chromosome name");
        }
        const std::string_view name(reinterpret_cast<const char*>(ptr), header.name_len);
        ptr += header.name_len;

        ChromosomeView chrom;
        chrom.name = name;
        chrom.length = header.length;
        chrom.text_length = header.text_length;
        chrom.primary_index = header.primary_index;
        chrom.occ_sample = occ_sample_;
        chrom.sa_sample = sa_sample_;
        std::copy(std::begin(header.c_array), std::end(header.c_array), chrom.c_array.begin());

        chrom.packed_bwt = mapping_ + header.packed_bwt_offset;
        chrom.packed_bwt_bytes = static_cast<std::size_t>(header.packed_bwt_bytes);
        chrom.sampled_occ = reinterpret_cast<const uint32_t*>(mapping_ + header.occ_offset);
        chrom.sampled_occ_entries = static_cast<std::size_t>(header.occ_entries);
        chrom.sampled_sa = reinterpret_cast<const uint32_t*>(mapping_ + header.sa_offset);
        chrom.sampled_sa_entries = static_cast<std::size_t>(header.sa_entries);
        
        // Set has_genome based on header flags
        if (is_v2) {
            chrom.has_genome = (header.flags & kFlagHasGenome) != 0;
        } else {
            chrom.has_genome = true;  // V1 always has genome
        }
        
        chrom.packed_genome = mapping_ + header.genome_offset;
        chrom.packed_genome_bytes = static_cast<std::size_t>(header.genome_bytes);

        if (header.genome_offset + header.genome_bytes > file_size_ ||
            header.sa_offset + header.sa_entries * sizeof(uint32_t) > file_size_ ||
            header.occ_offset + header.occ_entries * sizeof(uint32_t) > file_size_ ||
            header.packed_bwt_offset + header.packed_bwt_bytes > file_size_) {
            throw std::runtime_error("corrupt FM-index payload offsets");
        }
        chromosomes_.push_back(chrom);
    }
}

bool FMIndexView::is_open() const {
    return mapping_ != nullptr;
}

uint32_t FMIndexView::occ_stride() const {
    return occ_sample_;
}

uint32_t FMIndexView::sa_stride() const {
    return sa_sample_;
}

std::size_t FMIndexView::chromosome_count() const {
    return chromosomes_.size();
}

const FMIndexView::ChromosomeView& FMIndexView::chromosome(std::size_t index) const {
    return chromosomes_.at(index);
}

const std::vector<FMIndexView::ChromosomeView>& FMIndexView::chromosomes() const {
    return chromosomes_;
}

bool FMIndexView::has_genome() const {
    return (flags_ & kFlagHasGenome) != 0;
}

uint64_t FMIndexView::ChromosomeView::n() const {
    return text_length;
}

uint8_t FMIndexView::ChromosomeView::char_rank(uint64_t row) const {
    if (row == primary_index) {
        return kSentinelRank;
    }
    return get_packed_rank(packed_bwt, row);
}

uint32_t FMIndexView::ChromosomeView::occ(uint8_t rank, uint64_t pos) const {
    if (pos > text_length) {
        pos = text_length;
    }
    if (rank == kSentinelRank) {
        return static_cast<uint32_t>(pos > primary_index ? 1U : 0U);
    }
    if (rank == kSeparatorRank || rank > kTRank) {
        return 0U;
    }

    const uint64_t block = pos / occ_sample;
    const uint64_t start = block * occ_sample;
    uint32_t count = sampled_occ[static_cast<std::size_t>(block * 4ULL + dna_occ_index(rank))];
    count += simd::count_packed_range(static_cast<uint8_t>(rank - kARank), packed_bwt, start, pos, primary_index);
    return count;
}

uint64_t FMIndexView::ChromosomeView::lf(uint64_t row) const {
    const uint8_t rank = char_rank(row);
    return c_array[rank] + occ(rank, row);
}

uint64_t FMIndexView::ChromosomeView::locate(uint64_t row) const {
    uint64_t steps = 0;
    while (row % sa_sample != 0ULL) {
        row = lf(row);
        ++steps;
    }
    
    const std::size_t sa_idx = static_cast<std::size_t>(row / sa_sample);
    const uint32_t sa_value = sampled_sa[sa_idx];
    
    return (static_cast<uint64_t>(sa_value) + steps) % text_length;
}

bool FMIndexView::ChromosomeView::same_chromosome_window(uint64_t text_pos, uint64_t span) const {
    if (span == 0ULL || text_pos >= length) {
        return false;
    }
    return text_pos + span <= length;
}

void FMIndexView::ChromosomeView::extract_reference(uint64_t text_pos, uint64_t ref_length, std::string& out) const {
    out.resize(static_cast<std::size_t>(ref_length));
    for (uint64_t i = 0; i < ref_length; ++i) {
        out[static_cast<std::size_t>(i)] = base_from_rank(get_packed_rank(packed_genome, text_pos + i));
    }
}

}  // namespace mapper_memory

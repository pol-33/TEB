#include "index.hpp"

#include <algorithm>
#include <array>
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

namespace mapper_speed {

namespace {

constexpr char kMagic[8] = {'M', 'S', 'P', 'D', 'I', 'D', 'X', '1'};

const auto& packed_byte_to_bases() {
    static const std::array<uint32_t, 256> table = []() {
        std::array<uint32_t, 256> values{};
        for (uint32_t byte = 0; byte < 256u; ++byte) {
            char chars[4] = {
                kCodeToBase[byte & 0x3u],
                kCodeToBase[(byte >> 2u) & 0x3u],
                kCodeToBase[(byte >> 4u) & 0x3u],
                kCodeToBase[(byte >> 6u) & 0x3u]
            };
            uint32_t word = 0;
            std::memcpy(&word, chars, sizeof(word));
            values[byte] = word;
        }
        return values;
    }();
    return table;
}

template <typename T>
void write_value(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) {
        throw std::runtime_error("failed to write index file");
    }
}

void write_padding(std::ofstream& out, uint64_t target_offset) {
    const uint64_t current = static_cast<uint64_t>(out.tellp());
    if (current >= target_offset) {
        return;
    }
    std::vector<char> padding(static_cast<std::size_t>(target_offset - current), 0);
    out.write(padding.data(), static_cast<std::streamsize>(padding.size()));
    if (!out) {
        throw std::runtime_error("failed to write index padding");
    }
}

}  // namespace

std::size_t write_index(const std::string& path,
                        const ReferenceData& reference,
                        const std::vector<OffsetPageMeta>& page_meta,
                        const std::vector<uint8_t>& page_data,
                        const std::vector<uint32_t>& positions) {
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("failed to open index output: " + path);
    }

    IndexHeader header{};
    std::memcpy(header.magic, kMagic, sizeof(kMagic));
    header.version = kIndexVersion;
    header.checksum = reference.checksum;
    header.genome_length = reference.genome_length;
    header.chromosome_count = static_cast<uint32_t>(reference.chromosomes.size());
    header.packed_reference_bytes = reference.packed_bases.size();
    header.n_mask_bytes = reference.n_mask_words.size() * sizeof(uint64_t);
    header.offset_page_count = page_meta.size();
    header.positions_count = positions.size();

    uint64_t offset = sizeof(IndexHeader);
    header.chromosome_table_offset = offset;
    for (const auto& chrom : reference.chromosomes) {
        offset += sizeof(StoredChromosome) + chrom.name.size();
    }
    offset = align_up(offset, 64u);

    header.packed_reference_offset = offset;
    offset += header.packed_reference_bytes;
    offset = align_up(offset, 64u);

    header.n_mask_offset = offset;
    offset += header.n_mask_bytes;
    offset = align_up(offset, 64u);

    header.offset_page_meta_offset = offset;
    offset += static_cast<uint64_t>(page_meta.size()) * sizeof(OffsetPageMeta);
    offset = align_up(offset, 64u);

    const uint64_t page_data_offset = offset;
    offset += page_data.size();
    offset = align_up(offset, 64u);

    header.positions_offset = offset;
    offset += static_cast<uint64_t>(positions.size()) * sizeof(uint32_t);

    write_value(out, header);

    for (const auto& chrom : reference.chromosomes) {
        StoredChromosome stored{};
        stored.name_len = static_cast<uint32_t>(chrom.name.size());
        stored.start = chrom.start;
        stored.length = chrom.length;
        write_value(out, stored);
        out.write(chrom.name.data(), static_cast<std::streamsize>(chrom.name.size()));
        if (!out) {
            throw std::runtime_error("failed to write chromosome name");
        }
    }

    write_padding(out, header.packed_reference_offset);

    out.write(reinterpret_cast<const char*>(reference.packed_bases.data()),
              static_cast<std::streamsize>(reference.packed_bases.size()));
    write_padding(out, header.n_mask_offset);

    out.write(reinterpret_cast<const char*>(reference.n_mask_words.data()),
              static_cast<std::streamsize>(reference.n_mask_words.size() * sizeof(uint64_t)));
    write_padding(out, header.offset_page_meta_offset);

    std::vector<OffsetPageMeta> adjusted = page_meta;
    for (auto& meta : adjusted) {
        if ((meta.flags & (kOffsetPageDense | kOffsetPageSparse)) != 0u) {
            meta.data_offset += page_data_offset;
        }
    }
    out.write(reinterpret_cast<const char*>(adjusted.data()),
              static_cast<std::streamsize>(adjusted.size() * sizeof(OffsetPageMeta)));
    write_padding(out, page_data_offset);

    out.write(reinterpret_cast<const char*>(page_data.data()), static_cast<std::streamsize>(page_data.size()));
    write_padding(out, header.positions_offset);

    out.write(reinterpret_cast<const char*>(positions.data()),
              static_cast<std::streamsize>(positions.size() * sizeof(uint32_t)));
    if (!out) {
        throw std::runtime_error("failed to write index payload");
    }
    return static_cast<std::size_t>(out.tellp());
}

IndexView::IndexView() = default;

IndexView::IndexView(const std::string& path) {
    open(path);
}

IndexView::~IndexView() {
    close();
}

void IndexView::open(const std::string& path) {
    close();

    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("failed to open index: " + path);
    }

    struct stat st {};
    if (fstat(fd_, &st) != 0) {
        const int err = errno;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "failed to stat index");
    }
    mapping_size_ = static_cast<std::size_t>(st.st_size);
    mapping_ = static_cast<uint8_t*>(mmap(nullptr, mapping_size_, PROT_READ, MAP_PRIVATE, fd_, 0));
    if (mapping_ == MAP_FAILED) {
        const int err = errno;
        mapping_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "failed to mmap index");
    }

    std::memcpy(&header_, mapping_, sizeof(IndexHeader));
    if (std::memcmp(header_.magic, kMagic, sizeof(kMagic)) != 0) {
        throw std::runtime_error("invalid mapper-speed index magic");
    }
    if (header_.version != kIndexVersion || header_.seed_length != kSeedLength ||
        header_.page_shift != kOffsetPageShift || header_.flags == 0u) {
        throw std::runtime_error("unsupported mapper-speed index parameters");
    }

    const uint8_t* ptr = mapping_ + header_.chromosome_table_offset;
    chromosomes_.clear();
    chromosomes_.reserve(header_.chromosome_count);
    for (uint32_t i = 0; i < header_.chromosome_count; ++i) {
        StoredChromosome stored{};
        std::memcpy(&stored, ptr, sizeof(StoredChromosome));
        ptr += sizeof(StoredChromosome);
        ChromosomeRecord chrom;
        chrom.start = stored.start;
        chrom.length = stored.length;
        chrom.name.assign(reinterpret_cast<const char*>(ptr), stored.name_len);
        ptr += stored.name_len;
        chromosomes_.push_back(std::move(chrom));
    }

    packed_reference_ = mapping_ + header_.packed_reference_offset;
    n_mask_ = reinterpret_cast<const uint64_t*>(mapping_ + header_.n_mask_offset);
    page_meta_ = reinterpret_cast<const OffsetPageMeta*>(mapping_ + header_.offset_page_meta_offset);
    positions_ = reinterpret_cast<const uint32_t*>(mapping_ + header_.positions_offset);

    const std::size_t bin_count =
        (static_cast<std::size_t>(header_.genome_length) + kChromLookupBinSize - 1u) >> kChromLookupShift;
    chrom_lookup_bins_.assign(bin_count, 0);
    std::size_t chrom_index = 0;
    for (std::size_t bin = 0; bin < bin_count; ++bin) {
        const uint32_t pos = static_cast<uint32_t>(bin << kChromLookupShift);
        while (chrom_index + 1u < chromosomes_.size() && chromosomes_[chrom_index + 1u].start <= pos) {
            ++chrom_index;
        }
        chrom_lookup_bins_[bin] = static_cast<uint16_t>(chrom_index);
    }
}

void IndexView::close() {
    if (mapping_ != nullptr) {
        munmap(mapping_, mapping_size_);
        mapping_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
    mapping_size_ = 0;
    chromosomes_.clear();
    chrom_lookup_bins_.clear();
    packed_reference_ = nullptr;
    n_mask_ = nullptr;
    page_meta_ = nullptr;
    positions_ = nullptr;
}

bool IndexView::is_open() const {
    return mapping_ != nullptr;
}

uint32_t IndexView::genome_length() const {
    return header_.genome_length;
}

uint32_t IndexView::chromosome_count() const {
    return static_cast<uint32_t>(chromosomes_.size());
}

const ChromosomeRecord& IndexView::chromosome(std::size_t index) const {
    return chromosomes_.at(index);
}

uint32_t IndexView::positions_count() const {
    return static_cast<uint32_t>(header_.positions_count);
}

uint64_t IndexView::checksum() const {
    return header_.checksum;
}

uint32_t IndexView::index_stride() const {
    return header_.flags;
}

const OffsetPageMeta& IndexView::page_meta_for_entry(uint64_t entry) const {
    return page_meta_[entry >> kOffsetPageShift];
}

const uint32_t* IndexView::dense_page_values(const OffsetPageMeta& meta) const {
    return reinterpret_cast<const uint32_t*>(mapping_ + meta.data_offset);
}

const OffsetTransition* IndexView::sparse_page_records(const OffsetPageMeta& meta) const {
    return reinterpret_cast<const OffsetTransition*>(mapping_ + meta.data_offset);
}

uint32_t IndexView::offset_at(uint64_t key) const {
    const OffsetPageMeta& meta = page_meta_for_entry(key);
    if ((meta.flags & (kOffsetPageDense | kOffsetPageSparse)) == 0u || meta.record_count == 0u) {
        return meta.base_value;
    }
    const uint32_t local = static_cast<uint32_t>(key & (kOffsetPageSize - 1u));
    if ((meta.flags & kOffsetPageDense) != 0u) {
        const uint32_t* values = dense_page_values(meta);
        return values[local];
    }

    const OffsetTransition* records = sparse_page_records(meta);
    uint32_t lo = 0;
    uint32_t hi = meta.record_count;
    while (lo < hi) {
        const uint32_t mid = lo + ((hi - lo) >> 1u);
        if (records[mid].local_entry <= local) {
            lo = mid + 1u;
        } else {
            hi = mid;
        }
    }
    return (lo == 0u) ? meta.base_value : records[lo - 1u].value_after;
}

uint32_t IndexView::occurrence_count(uint32_t key) const {
    return offset_at(static_cast<uint64_t>(key) + 1u) - offset_at(key);
}

std::pair<const uint32_t*, const uint32_t*> IndexView::positions_for(uint32_t key) const {
    const uint32_t begin = offset_at(key);
    const uint32_t end = offset_at(static_cast<uint64_t>(key) + 1u);
    return {positions_ + begin, positions_ + end};
}

bool IndexView::has_n(uint32_t global_pos) const {
    return bitset_get(n_mask_, global_pos);
}

char IndexView::base_at(uint32_t global_pos) const {
    if (has_n(global_pos)) {
        return 'N';
    }
    return kCodeToBase[packed_get(packed_reference_, global_pos)];
}

void IndexView::extract_sequence(uint32_t global_pos, uint32_t length, std::string& out) const {
    out.resize(length);
    if (length == 0u) {
        return;
    }

    char* dst = out.data();
    uint32_t pos = global_pos;
    uint32_t written = 0;

    const uint32_t prefix = std::min<uint32_t>(length, static_cast<uint32_t>((4u - (pos & 3u)) & 3u));
    for (; written < prefix; ++written, ++pos) {
        dst[written] = kCodeToBase[packed_get(packed_reference_, pos)];
    }

    const auto& decode = packed_byte_to_bases();
    const uint32_t full_bytes = (length - written) >> 2u;
    std::size_t byte_index = pos >> 2u;
    for (uint32_t i = 0; i < full_bytes; ++i) {
        const uint32_t packed_chars = decode[packed_reference_[byte_index + i]];
        std::memcpy(dst + written + (i << 2u), &packed_chars, sizeof(uint32_t));
    }
    written += full_bytes << 2u;
    pos += full_bytes << 2u;

    for (; written < length; ++written, ++pos) {
        dst[written] = kCodeToBase[packed_get(packed_reference_, pos)];
    }

    const uint32_t end_pos = global_pos + length;
    std::size_t word_index = global_pos >> 6u;
    const std::size_t end_word = (end_pos - 1u) >> 6u;
    while (word_index <= end_word) {
        uint64_t mask = n_mask_[word_index];
        const uint32_t word_base = static_cast<uint32_t>(word_index << 6u);
        if (word_base < global_pos) {
            mask &= (~uint64_t{0}) << (global_pos - word_base);
        }
        if (word_base + 64u > end_pos) {
            mask &= (~uint64_t{0}) >> ((word_base + 64u) - end_pos);
        }
        while (mask != 0u) {
            const uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(mask));
            dst[word_base + bit - global_pos] = 'N';
            mask &= (mask - 1u);
        }
        ++word_index;
    }
}

uint32_t IndexView::count_mismatches_packed(uint32_t global_pos,
                                            const uint8_t* packed_read,
                                            uint32_t length,
                                            uint32_t limit) const {
    constexpr uint64_t kLaneMask = UINT64_C(0x5555555555555555);
    uint32_t mismatches = 0;
    uint32_t offset = 0;

    auto load_ref_codes = [&](uint32_t pos) -> uint64_t {
        const std::size_t byte_index = pos >> 2u;
        const unsigned bit_shift = (pos & 3u) * 2u;
        uint64_t lo = 0;
        uint64_t hi = 0;
        const std::size_t available = header_.packed_reference_bytes - byte_index;
        const std::size_t lo_bytes = std::min<std::size_t>(sizeof(lo), available);
        if (lo_bytes > 0u) {
            std::memcpy(&lo, packed_reference_ + byte_index, lo_bytes);
        }
        if (bit_shift == 0u) {
            return lo;
        }
        if (available > sizeof(lo)) {
            const std::size_t hi_bytes = std::min<std::size_t>(sizeof(hi), available - sizeof(lo));
            std::memcpy(&hi, packed_reference_ + byte_index + sizeof(lo), hi_bytes);
        }
        return (lo >> bit_shift) | (hi << (64u - bit_shift));
    };

    auto load_ref_n_mask = [&](uint32_t pos) -> uint64_t {
        const std::size_t word_index = pos >> 6u;
        const unsigned bit_shift = pos & 63u;
        uint64_t lo = n_mask_[word_index];
        if (bit_shift == 0u) {
            return lo;
        }
        uint64_t hi = 0;
        const std::size_t next_word = word_index + 1u;
        const std::size_t n_mask_words = header_.n_mask_bytes / sizeof(uint64_t);
        if (next_word < n_mask_words) {
            hi = n_mask_[next_word];
        }
        return (lo >> bit_shift) | (hi << (64u - bit_shift));
    };

    for (; offset + 32u <= length; offset += 32u) {
        uint64_t read_chunk = 0;
        std::memcpy(&read_chunk, packed_read + (offset >> 2u), sizeof(read_chunk));
        const uint64_t ref_chunk = load_ref_codes(global_pos + offset);
        const uint32_t n_chunk = static_cast<uint32_t>(load_ref_n_mask(global_pos + offset));
        if (n_chunk == 0u) {
            const uint64_t diff = read_chunk ^ ref_chunk;
            mismatches += static_cast<uint32_t>(__builtin_popcountll((diff | (diff >> 1u)) & kLaneMask));
        } else {
            for (uint32_t i = 0; i < 32u; ++i) {
                if (((n_chunk >> i) & 1u) != 0u) {
                    ++mismatches;
                    continue;
                }
                const uint8_t ref_code = static_cast<uint8_t>((ref_chunk >> (i * 2u)) & 0x3u);
                const uint8_t read_code = static_cast<uint8_t>((read_chunk >> (i * 2u)) & 0x3u);
                mismatches += static_cast<uint32_t>(ref_code != read_code);
            }
        }
        if (mismatches > limit) {
            return mismatches;
        }
    }

    for (; offset < length; ++offset) {
        const uint32_t pos = global_pos + offset;
        if (has_n(pos) ||
            packed_get(packed_read, offset) != packed_get(packed_reference_, pos)) {
            ++mismatches;
            if (mismatches > limit) {
                return mismatches;
            }
        }
    }
    return mismatches;
}

std::size_t IndexView::chromosome_for_position(uint32_t global_pos) const {
    if (chrom_lookup_bins_.empty()) {
        return 0u;
    }
    std::size_t index =
        chrom_lookup_bins_[std::min<std::size_t>(global_pos >> kChromLookupShift, chrom_lookup_bins_.size() - 1u)];
    while (index + 1u < chromosomes_.size() && chromosomes_[index + 1u].start <= global_pos) {
        ++index;
    }
    while (index > 0u && chromosomes_[index].start > global_pos) {
        --index;
    }
    return index;
}

bool IndexView::stays_within_chromosome(uint32_t global_pos, uint32_t ref_length) const {
    const std::size_t chrom_index = chromosome_for_position(global_pos);
    const ChromosomeRecord& chrom = chromosomes_[chrom_index];
    const uint64_t local = static_cast<uint64_t>(global_pos) - chrom.start;
    return local + ref_length <= chrom.length;
}

}  // namespace mapper_speed

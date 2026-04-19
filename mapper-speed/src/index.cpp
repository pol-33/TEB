#include "index.hpp"

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

namespace mapper_speed {

namespace {

constexpr char kMagic[8] = {'M', 'S', 'P', 'D', 'I', 'D', 'X', '1'};

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
    header.flags = kIndexStride;
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
        header_.page_shift != kOffsetPageShift || header_.flags != kIndexStride) {
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
    for (uint32_t i = 0; i < length; ++i) {
        out[i] = base_at(global_pos + i);
    }
}

std::size_t IndexView::chromosome_for_position(uint32_t global_pos) const {
    auto it = std::upper_bound(chromosomes_.begin(), chromosomes_.end(), global_pos,
        [](uint32_t pos, const ChromosomeRecord& chrom) {
            return pos < chrom.start;
        });
    if (it == chromosomes_.begin()) {
        return 0;
    }
    --it;
    return static_cast<std::size_t>(it - chromosomes_.begin());
}

bool IndexView::stays_within_chromosome(uint32_t global_pos, uint32_t ref_length) const {
    const std::size_t chrom_index = chromosome_for_position(global_pos);
    const ChromosomeRecord& chrom = chromosomes_[chrom_index];
    const uint64_t local = static_cast<uint64_t>(global_pos) - chrom.start;
    return local + ref_length <= chrom.length;
}

}  // namespace mapper_speed

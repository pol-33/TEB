#ifndef MAPPER_SPEED_INDEX_HPP
#define MAPPER_SPEED_INDEX_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "common.hpp"
#include "reference.hpp"

namespace mapper_speed {

struct OffsetPageMeta {
    uint64_t data_offset = 0;
    uint32_t base_value = 0;
    uint16_t value_count = 0;
    uint8_t flags = 0;
    uint8_t reserved = 0;
};

enum OffsetPageFlags : uint8_t {
    kOffsetPageDense = 1
};

struct IndexHeader {
    char magic[8];
    uint32_t version = 1;
    uint32_t seed_length = kSeedLength;
    uint32_t page_shift = kOffsetPageShift;
    uint32_t flags = 0;
    uint64_t checksum = 0;
    uint32_t genome_length = 0;
    uint32_t chromosome_count = 0;
    uint64_t chromosome_table_offset = 0;
    uint64_t packed_reference_offset = 0;
    uint64_t packed_reference_bytes = 0;
    uint64_t n_mask_offset = 0;
    uint64_t n_mask_bytes = 0;
    uint64_t offset_page_meta_offset = 0;
    uint64_t offset_page_count = 0;
    uint64_t positions_offset = 0;
    uint64_t positions_count = 0;
    uint64_t high_freq_bitmap_offset = 0;
    uint64_t high_freq_bitmap_bytes = 0;
};

struct StoredChromosome {
    uint32_t name_len = 0;
    uint32_t start = 0;
    uint32_t length = 0;
};

struct DenseBuildOutput {
    std::vector<OffsetPageMeta> page_meta;
    std::vector<uint8_t> dense_pages;
    std::vector<uint32_t> positions;
};

std::size_t write_index(const std::string& path,
                        const ReferenceData& reference,
                        const std::vector<OffsetPageMeta>& page_meta,
                        const std::vector<uint8_t>& dense_pages,
                        const std::vector<uint32_t>& positions);

class IndexView {
public:
    IndexView();
    explicit IndexView(const std::string& path);
    ~IndexView();

    IndexView(const IndexView&) = delete;
    IndexView& operator=(const IndexView&) = delete;

    void open(const std::string& path);
    void close();
    bool is_open() const;

    uint32_t genome_length() const;
    uint32_t chromosome_count() const;
    const ChromosomeRecord& chromosome(std::size_t index) const;
    uint32_t positions_count() const;
    uint64_t checksum() const;

    uint32_t offset_at(uint64_t key) const;
    uint32_t occurrence_count(uint32_t key) const;
    std::pair<const uint32_t*, const uint32_t*> positions_for(uint32_t key) const;
    bool has_n(uint32_t global_pos) const;
    char base_at(uint32_t global_pos) const;
    void extract_sequence(uint32_t global_pos, uint32_t length, std::string& out) const;
    std::size_t chromosome_for_position(uint32_t global_pos) const;
    bool stays_within_chromosome(uint32_t global_pos, uint32_t ref_length) const;

private:
    const uint32_t* dense_page_values(const OffsetPageMeta& meta) const;
    const OffsetPageMeta& page_meta_for_entry(uint64_t entry) const;

    int fd_ = -1;
    std::size_t mapping_size_ = 0;
    uint8_t* mapping_ = nullptr;
    IndexHeader header_{};
    std::vector<ChromosomeRecord> chromosomes_;
    const uint8_t* packed_reference_ = nullptr;
    const uint64_t* n_mask_ = nullptr;
    const OffsetPageMeta* page_meta_ = nullptr;
    const uint32_t* positions_ = nullptr;
};

}  // namespace mapper_speed

#endif

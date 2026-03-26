#ifndef MAPPER_MEMORY_SEED_INDEX_HPP
#define MAPPER_MEMORY_SEED_INDEX_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "fasta_reader.hpp"

namespace mapper_memory {

constexpr uint32_t kSeedLength = 13;
constexpr uint64_t kSeedBucketCount = 1ULL << (2U * kSeedLength);

struct IndexBuildStats {
    uint64_t genome_length = 0;
    uint64_t indexed_positions = 0;
    uint32_t num_chromosomes = 0;
    uint64_t file_size = 0;
};

uint32_t encode_seed(std::string_view seed);

IndexBuildStats build_seed_index(const std::string& fasta_path, const std::string& index_path);

class SeedIndexView {
public:
    struct ChromView {
        std::string_view name;
        uint64_t offset = 0;
        uint64_t length = 0;
    };

    SeedIndexView();
    explicit SeedIndexView(const std::string& path);
    ~SeedIndexView();

    SeedIndexView(const SeedIndexView&) = delete;
    SeedIndexView& operator=(const SeedIndexView&) = delete;
    SeedIndexView(SeedIndexView&& other) noexcept;
    SeedIndexView& operator=(SeedIndexView&& other) noexcept;

    void open(const std::string& path);
    void close();
    bool is_open() const;

    uint32_t seed_length() const;
    uint64_t genome_length() const;
    uint64_t indexed_positions() const;
    const std::vector<ChromView>& chromosomes() const;

    uint32_t bucket_begin(uint32_t code) const;
    uint32_t bucket_end(uint32_t code) const;
    uint32_t position_at(uint64_t idx) const;

    uint64_t chromosome_index(uint64_t genome_pos) const;
    bool same_chromosome_window(uint64_t genome_pos, uint64_t span) const;
    std::string_view chromosome_name(uint64_t chrom_index) const;
    uint64_t chromosome_offset(uint64_t chrom_index) const;
    uint64_t chromosome_length(uint64_t chrom_index) const;

    char reference_base(uint64_t genome_pos) const;
    void extract_reference(uint64_t genome_pos, uint64_t length, std::string& out) const;

private:
    void parse();
    void move_from(SeedIndexView&& other) noexcept;

    int fd_;
    const uint8_t* mapping_;
    std::size_t file_size_;

    uint32_t seed_length_;
    uint64_t genome_length_;
    uint64_t indexed_positions_;
    const uint8_t* packed_genome_;
    const uint32_t* offsets_;
    const uint32_t* positions_;
    std::vector<ChromView> chromosomes_;
};

}  // namespace mapper_memory

#endif

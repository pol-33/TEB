#ifndef MAPPER_MEMORY_FM_INDEX_HPP
#define MAPPER_MEMORY_FM_INDEX_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "bwt.hpp"
#include "fasta_reader.hpp"

namespace mapper_memory {

class OwnedFMIndex {
public:
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    uint32_t occ_sample = 128;
    uint32_t sa_sample = 32;
    std::vector<ChromInfo> chromosomes;
    std::array<uint32_t, 5> c_array{};
    std::vector<uint8_t> packed_bwt;
    std::vector<uint32_t> sampled_occ;
    std::vector<uint32_t> sampled_sa;
    std::vector<uint8_t> packed_genome;

    std::size_t write(const std::string& path) const;
};

OwnedFMIndex build_fm_index(const std::string& genome,
                            const std::vector<ChromInfo>& chromosomes,
                            const BWTData& bwt,
                            uint32_t occ_sample = 128,
                            uint32_t sa_sample = 32);

class FMIndexView {
public:
    struct ChromView {
        std::string_view name;
        uint64_t offset = 0;
        uint64_t length = 0;
    };

    FMIndexView();
    explicit FMIndexView(const std::string& path);
    ~FMIndexView();

    FMIndexView(const FMIndexView&) = delete;
    FMIndexView& operator=(const FMIndexView&) = delete;
    FMIndexView(FMIndexView&& other) noexcept;
    FMIndexView& operator=(FMIndexView&& other) noexcept;

    void open(const std::string& path);
    void close();

    bool is_open() const;
    uint64_t n() const;
    uint64_t genome_length() const;
    uint64_t primary() const;
    uint32_t occ_stride() const;
    uint32_t sa_stride() const;
    const std::array<uint32_t, 5>& c_array() const;
    const std::vector<ChromView>& chromosomes() const;

    uint8_t char_code(uint64_t row) const;
    uint32_t occ(uint8_t code, uint64_t pos) const;
    uint64_t lf(uint64_t row) const;
    uint64_t sa_value(uint64_t row) const;
    uint64_t chromosome_index(uint64_t genome_pos) const;
    bool same_chromosome_window(uint64_t genome_pos, uint64_t span) const;
    std::string_view chromosome_name(uint64_t chrom_index) const;
    uint64_t chromosome_offset(uint64_t chrom_index) const;
    uint64_t chromosome_length(uint64_t chrom_index) const;
    char reference_base(uint64_t genome_pos) const;
    void extract_reference(uint64_t genome_pos, uint64_t length, std::string& out) const;

private:
    void parse();
    void move_from(FMIndexView&& other) noexcept;

    int fd_;
    const uint8_t* mapping_;
    std::size_t file_size_;

    uint64_t text_length_;
    uint64_t primary_index_;
    uint32_t occ_sample_;
    uint32_t sa_sample_;
    std::array<uint32_t, 5> c_array_values_;
    const uint8_t* packed_bwt_;
    const uint32_t* sampled_occ_;
    const uint32_t* sampled_sa_;
    const uint8_t* packed_genome_;
    std::vector<ChromView> chromosomes_;
};

}  // namespace mapper_memory

#endif

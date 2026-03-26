#ifndef MAPPER_MEMORY_FM_INDEX_HPP
#define MAPPER_MEMORY_FM_INDEX_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "bwt.hpp"

namespace mapper_memory {

struct OwnedChromosomeIndex {
    std::string name;
    uint64_t length = 0;
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    uint32_t occ_sample = 256;
    uint32_t sa_sample = 32;
    std::array<uint64_t, kAlphabetSize> c_array{};
    std::vector<uint8_t> packed_bwt;
    std::vector<uint32_t> sampled_occ;  // row-major [block][A,C,G,T]
    std::vector<uint32_t> sampled_sa;
    std::vector<uint8_t> packed_genome;
};

OwnedChromosomeIndex build_chromosome_index(const std::string& name,
                                            const std::string& sequence,
                                            const BWTData& bwt,
                                            const std::vector<uint32_t>& suffix_array,
                                            uint32_t occ_sample = 256,
                                            uint32_t sa_sample = 32);

class FMIndexView {
public:
    struct ChromosomeView {
        std::string_view name;
        uint64_t length = 0;
        uint64_t text_length = 0;
        uint64_t primary_index = 0;
        uint32_t occ_sample = 0;
        uint32_t sa_sample = 0;
        std::array<uint64_t, kAlphabetSize> c_array{};
        const uint8_t* packed_bwt = nullptr;
        std::size_t packed_bwt_bytes = 0;
        const uint32_t* sampled_occ = nullptr;
        std::size_t sampled_occ_entries = 0;
        const uint32_t* sampled_sa = nullptr;
        std::size_t sampled_sa_entries = 0;
        const uint8_t* packed_genome = nullptr;
        std::size_t packed_genome_bytes = 0;

        uint64_t n() const;
        uint8_t char_rank(uint64_t row) const;
        uint32_t occ(uint8_t rank, uint64_t pos) const;
        uint64_t lf(uint64_t row) const;
        uint64_t locate(uint64_t row) const;
        bool same_chromosome_window(uint64_t text_pos, uint64_t span) const;
        void extract_reference(uint64_t text_pos, uint64_t length, std::string& out) const;
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

    uint32_t occ_stride() const;
    uint32_t sa_stride() const;
    std::size_t chromosome_count() const;
    const ChromosomeView& chromosome(std::size_t index) const;
    const std::vector<ChromosomeView>& chromosomes() const;

    static std::size_t write(const std::string& path,
                             const std::vector<OwnedChromosomeIndex>& chromosomes,
                             uint32_t occ_sample = 256,
                             uint32_t sa_sample = 32);

private:
    void parse();
    void move_from(FMIndexView&& other) noexcept;

    int fd_;
    const uint8_t* mapping_;
    std::size_t file_size_;
    uint32_t occ_sample_;
    uint32_t sa_sample_;
    std::vector<ChromosomeView> chromosomes_;
};

}  // namespace mapper_memory

#endif

#ifndef MAPPER_MEMORY_BWT_HPP
#define MAPPER_MEMORY_BWT_HPP

#include <array>
#include <cstdint>
#include <vector>

#include "nucleotide.hpp"

namespace mapper_memory {

struct BWTData {
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    std::array<uint64_t, kAlphabetSize> counts{};
    std::vector<uint8_t> packed_bwt;
    std::vector<uint32_t> separator_rows;
};

std::vector<uint32_t> build_suffix_array_sais(const std::vector<uint8_t>& text);
BWTData build_bwt(const std::vector<uint8_t>& text, const std::vector<uint32_t>& suffix_array);

}  // namespace mapper_memory

#endif

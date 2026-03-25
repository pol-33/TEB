#ifndef MAPPER_MEMORY_BWT_HPP
#define MAPPER_MEMORY_BWT_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace mapper_memory {

struct BWTData {
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    std::vector<uint8_t> packed_bwt;
    std::vector<uint32_t> suffix_array;
};

std::vector<uint32_t> build_suffix_array_sais(const std::string& genome);
BWTData build_bwt(const std::string& genome, const std::vector<uint32_t>& suffix_array);

}  // namespace mapper_memory

#endif

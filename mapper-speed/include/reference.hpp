#ifndef MAPPER_SPEED_REFERENCE_HPP
#define MAPPER_SPEED_REFERENCE_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace mapper_speed {

struct ChromosomeRecord {
    std::string name;
    uint32_t start = 0;
    uint32_t length = 0;
};

struct ReferenceData {
    std::vector<ChromosomeRecord> chromosomes;
    std::vector<uint8_t> packed_bases;
    std::vector<uint64_t> n_mask_words;
    uint32_t genome_length = 0;
    uint64_t checksum = 0;
};

ReferenceData load_reference(const std::string& path);

}  // namespace mapper_speed

#endif

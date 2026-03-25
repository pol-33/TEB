#ifndef MAPPER_MEMORY_FASTA_READER_HPP
#define MAPPER_MEMORY_FASTA_READER_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace mapper_memory {

struct ChromInfo {
    std::string name;
    uint64_t offset = 0;
    uint64_t length = 0;
};

struct FastaData {
    std::string genome;
    std::vector<ChromInfo> chromosomes;
};

FastaData load_fasta(const std::string& path);

}  // namespace mapper_memory

#endif

#ifndef MAPPER_MEMORY_FASTA_READER_HPP
#define MAPPER_MEMORY_FASTA_READER_HPP

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

namespace mapper_memory {

struct ChromInfo {
    std::string name;
    uint64_t offset = 0;
    uint64_t length = 0;
};

struct FastaData {
    uint64_t genome_length = 0;
    std::vector<uint8_t> packed_genome;
    std::vector<uint8_t> text;
    std::vector<ChromInfo> chromosomes;

    uint64_t text_length() const {
        return static_cast<uint64_t>(text.size());
    }
};

struct FastaChromosome {
    std::string name;
    uint64_t offset = 0;
    std::string sequence;
};

class FastaReader {
public:
    explicit FastaReader(const std::string& path);

    bool next(FastaChromosome& chrom);

private:
    std::ifstream in_;
    std::string pending_header_;
    uint64_t next_offset_ = 0;
};

FastaData load_fasta(const std::string& path);

}  // namespace mapper_memory

#endif

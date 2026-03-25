#ifndef MAPPER_MEMORY_FASTQ_READER_HPP
#define MAPPER_MEMORY_FASTQ_READER_HPP

#include <fstream>
#include <string>

namespace mapper_memory {

struct Read {
    std::string name;
    std::string seq;
    std::string qual;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path);

    bool next(Read& read);
    bool good() const;

private:
    std::ifstream in_;
};

}  // namespace mapper_memory

#endif

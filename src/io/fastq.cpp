#include "fastq.hpp"
#include <stdexcept>

FastqReader::FastqReader(const std::string&) {
    throw std::runtime_error("not implemented");
}

bool FastqReader::next(Read&) {
    throw std::runtime_error("not implemented");
}

#include "fastq_reader.hpp"

#include <stdexcept>

namespace mapper_memory {

FastqReader::FastqReader(const std::string& path) : in_(path) {
    if (!in_) {
        throw std::runtime_error("failed to open FASTQ file: " + path);
    }
}

bool FastqReader::next(Read& read) {
    std::string header;
    std::string plus;
    if (!std::getline(in_, header)) {
        return false;
    }
    if (!std::getline(in_, read.seq) || !std::getline(in_, plus) || !std::getline(in_, read.qual)) {
        throw std::runtime_error("truncated FASTQ record");
    }
    if (!header.empty() && header.back() == '\r') {
        header.pop_back();
    }
    if (!read.seq.empty() && read.seq.back() == '\r') {
        read.seq.pop_back();
    }
    if (!plus.empty() && plus.back() == '\r') {
        plus.pop_back();
    }
    if (!read.qual.empty() && read.qual.back() == '\r') {
        read.qual.pop_back();
    }
    if (header.empty() || header.front() != '@') {
        throw std::runtime_error("invalid FASTQ header");
    }
    if (plus.empty() || plus.front() != '+') {
        throw std::runtime_error("invalid FASTQ separator");
    }
    if (read.seq.size() != read.qual.size()) {
        throw std::runtime_error("FASTQ sequence/quality length mismatch");
    }
    read.name = header.substr(1);
    return true;
}

bool FastqReader::good() const {
    return in_.good();
}

}  // namespace mapper_memory

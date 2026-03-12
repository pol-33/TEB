#include "kmer_index.hpp"
#include <stdexcept>

KmerIndex::KmerIndex(size_t k) : k_(k), table_() {}

void KmerIndex::build(const std::string&) {
    throw std::runtime_error("not implemented");
}

void KmerIndex::save(const std::string&) const {
    throw std::runtime_error("not implemented");
}

void KmerIndex::load(const std::string&, LoadMode) {
    throw std::runtime_error("not implemented");
}

std::vector<uint32_t> KmerIndex::query(const std::string&) const {
    throw std::runtime_error("not implemented");
}

int KmerIndex::verify_match(uint32_t, const uint8_t*, uint32_t) const {
    throw std::runtime_error("not implemented");
}

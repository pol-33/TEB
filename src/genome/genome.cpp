#include "genome.hpp"
#include <istream>
#include <ostream>
#include <stdexcept>

void GenomeStorage::encode(const std::string&) {
    throw std::runtime_error("not implemented");
}

char GenomeStorage::get_base(uint32_t) const {
    throw std::runtime_error("not implemented");
}

std::string GenomeStorage::substr(uint32_t, uint32_t) const {
    throw std::runtime_error("not implemented");
}

bool GenomeStorage::is_n(uint32_t) const {
    throw std::runtime_error("not implemented");
}

uint32_t GenomeStorage::size() const {
    throw std::runtime_error("not implemented");
}

const uint8_t* GenomeStorage::packed_data() const {
    throw std::runtime_error("not implemented");
}

const uint64_t* GenomeStorage::n_mask_data() const {
    throw std::runtime_error("not implemented");
}

void GenomeStorage::write_binary(std::ostream&) const {
    throw std::runtime_error("not implemented");
}

void GenomeStorage::read_binary(std::istream&) {
    throw std::runtime_error("not implemented");
}

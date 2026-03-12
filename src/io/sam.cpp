#include "sam.hpp"
#include <stdexcept>

SamWriter::SamWriter(const std::string&) {
    throw std::runtime_error("not implemented");
}

void SamWriter::write_header(const std::string&, uint32_t) {
    throw std::runtime_error("not implemented");
}

void SamWriter::write_alignment(const std::string&, uint32_t,
                                const std::string&, int, bool) {
    throw std::runtime_error("not implemented");
}

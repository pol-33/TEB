#pragma once

#include <ostream>
#include <istream>
#include <stdexcept>
#include <vector>

namespace bio {

template<typename T>
void write_pod(std::ostream&, const T&) {
    throw std::runtime_error("not implemented");
}

template<typename T>
void read_pod(std::istream&, T&) {
    throw std::runtime_error("not implemented");
}

template<typename T>
void write_vec(std::ostream&, const std::vector<T>&) {
    throw std::runtime_error("not implemented");
}

template<typename T>
void read_vec(std::istream&, std::vector<T>&) {
    throw std::runtime_error("not implemented");
}

} // namespace bio

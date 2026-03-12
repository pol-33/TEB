#pragma once

#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace bio {

// Write a single trivially-copyable value as raw bytes.
template<typename T>
void write_pod(std::ostream& os, const T& val) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "bio::write_pod requires a trivially copyable type");
    os.write(reinterpret_cast<const char*>(&val), sizeof(T));
    if (!os) throw std::runtime_error("bio::write_pod: write failed");
}

// Read a single trivially-copyable value from raw bytes.
template<typename T>
void read_pod(std::istream& is, T& val) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "bio::read_pod requires a trivially copyable type");
    is.read(reinterpret_cast<char*>(&val), sizeof(T));
    if (!is) throw std::runtime_error("bio::read_pod: read failed");
}

// Write format: [count : uint64_t (8B)][element data : count * sizeof(T) B]
template<typename T>
void write_vec(std::ostream& os, const std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "bio::write_vec requires trivially copyable element type");
    const uint64_t n = static_cast<uint64_t>(v.size());
    write_pod(os, n);
    if (n > 0)
        os.write(reinterpret_cast<const char*>(v.data()),
                 static_cast<std::streamsize>(n * sizeof(T)));
    if (!os) throw std::runtime_error("bio::write_vec: write failed");
}

// Read format: [count : uint64_t (8B)][element data : count * sizeof(T) B]
template<typename T>
void read_vec(std::istream& is, std::vector<T>& v) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "bio::read_vec requires trivially copyable element type");
    uint64_t n = 0;
    read_pod(is, n);
    v.resize(static_cast<typename std::vector<T>::size_type>(n));
    if (n > 0)
        is.read(reinterpret_cast<char*>(v.data()),
                static_cast<std::streamsize>(n * sizeof(T)));
    if (!is) throw std::runtime_error("bio::read_vec: read failed");
}

} // namespace bio

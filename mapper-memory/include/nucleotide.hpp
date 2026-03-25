#ifndef MAPPER_MEMORY_NUCLEOTIDE_HPP
#define MAPPER_MEMORY_NUCLEOTIDE_HPP

#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace mapper_memory {

constexpr uint8_t kSentinelCode = 0;
constexpr uint8_t kACode = 1;
constexpr uint8_t kCCode = 2;
constexpr uint8_t kGCode = 3;
constexpr uint8_t kTCode = 4;

inline char normalize_base(char base) {
    switch (std::toupper(static_cast<unsigned char>(base))) {
        case 'A':
            return 'A';
        case 'C':
            return 'C';
        case 'G':
            return 'G';
        case 'T':
            return 'T';
        case 'N':
            return 'A';
        default:
            return 'A';
    }
}

inline uint8_t encode_base(char base) {
    switch (normalize_base(base)) {
        case 'A':
            return kACode;
        case 'C':
            return kCCode;
        case 'G':
            return kGCode;
        case 'T':
            return kTCode;
        default:
            return kACode;
    }
}

inline char decode_base(uint8_t code) {
    switch (code) {
        case kSentinelCode:
            return '$';
        case kACode:
            return 'A';
        case kCCode:
            return 'C';
        case kGCode:
            return 'G';
        case kTCode:
            return 'T';
        default:
            return 'A';
    }
}

inline std::string normalize_sequence(const std::string& sequence) {
    std::string out;
    out.reserve(sequence.size());
    for (char base : sequence) {
        out.push_back(normalize_base(base));
    }
    return out;
}

inline uint8_t packed_base_from_code(uint8_t code) {
    if (code < kACode || code > kTCode) {
        throw std::runtime_error("packed_base_from_code received a non-base code");
    }
    return static_cast<uint8_t>(code - 1);
}

inline uint8_t code_from_packed_base(uint8_t packed) {
    return static_cast<uint8_t>(packed + 1);
}

inline std::size_t packed_byte_count(uint64_t symbols) {
    return static_cast<std::size_t>((symbols + 3ULL) / 4ULL);
}

inline void set_packed_code(std::vector<uint8_t>& packed, uint64_t index, uint8_t code) {
    const std::size_t byte_index = static_cast<std::size_t>(index / 4ULL);
    const uint8_t shift = static_cast<uint8_t>((index % 4ULL) * 2ULL);
    const uint8_t value = packed_base_from_code(code);
    packed[byte_index] &= static_cast<uint8_t>(~(0x3u << shift));
    packed[byte_index] |= static_cast<uint8_t>(value << shift);
}

inline uint8_t get_packed_code(const uint8_t* packed, uint64_t index) {
    const std::size_t byte_index = static_cast<std::size_t>(index / 4ULL);
    const uint8_t shift = static_cast<uint8_t>((index % 4ULL) * 2ULL);
    const uint8_t packed_value = static_cast<uint8_t>((packed[byte_index] >> shift) & 0x3u);
    return code_from_packed_base(packed_value);
}

inline std::vector<uint8_t> pack_sequence(const std::string& sequence) {
    std::vector<uint8_t> packed(packed_byte_count(sequence.size()), 0);
    for (uint64_t i = 0; i < sequence.size(); ++i) {
        set_packed_code(packed, i, encode_base(sequence[static_cast<std::size_t>(i)]));
    }
    return packed;
}

}  // namespace mapper_memory

#endif

#ifndef MAPPER_SPEED_NUCLEOTIDE_HPP
#define MAPPER_SPEED_NUCLEOTIDE_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace mapper_speed {

inline constexpr uint8_t kBaseCodeInvalid = 0xFF;

inline constexpr uint8_t kAsciiToCode[32] = {
    kBaseCodeInvalid, 0, kBaseCodeInvalid, 1, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, 2,
    kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid,
    kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, 3, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid,
    kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid, kBaseCodeInvalid
};

inline constexpr char kCodeToBase[4] = {'A', 'C', 'G', 'T'};

inline uint8_t base_to_code(char base) {
    return kAsciiToCode[(static_cast<uint8_t>(base) & static_cast<uint8_t>(~0x20)) & 0x1F];
}

inline char normalize_base(char base) {
    const uint8_t code = base_to_code(base);
    return (code == kBaseCodeInvalid) ? 'N' : kCodeToBase[code];
}

inline char complement_base(char base) {
    switch (normalize_base(base)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'N';
    }
}

inline std::string normalize_sequence(const std::string& sequence) {
    std::string out(sequence.size(), 'N');
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        out[i] = normalize_base(sequence[i]);
    }
    return out;
}

inline std::string reverse_complement(const std::string& sequence) {
    std::string out(sequence.size(), 'N');
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        out[i] = complement_base(sequence[sequence.size() - 1 - i]);
    }
    return out;
}

inline std::size_t packed_base_bytes(std::size_t n) {
    return (n + 3u) / 4u;
}

inline void packed_set(std::vector<uint8_t>& packed, std::size_t index, uint8_t code) {
    const std::size_t byte_index = index >> 2u;
    if (packed.size() <= byte_index) {
        packed.resize(byte_index + 1u, 0);
    }
    const uint8_t shift = static_cast<uint8_t>((index & 3u) * 2u);
    packed[byte_index] &= static_cast<uint8_t>(~(0x3u << shift));
    packed[byte_index] |= static_cast<uint8_t>((code & 0x3u) << shift);
}

inline uint8_t packed_get(const uint8_t* packed, std::size_t index) {
    const std::size_t byte_index = index >> 2u;
    const uint8_t shift = static_cast<uint8_t>((index & 3u) * 2u);
    return static_cast<uint8_t>((packed[byte_index] >> shift) & 0x3u);
}

inline void bitset_set(std::vector<uint64_t>& words, std::size_t index) {
    const std::size_t word_index = index >> 6u;
    if (words.size() <= word_index) {
        words.resize(word_index + 1u, 0);
    }
    words[word_index] |= (uint64_t{1} << (index & 63u));
}

inline bool bitset_get(const uint64_t* words, std::size_t index) {
    return ((words[index >> 6u] >> (index & 63u)) & 1ULL) != 0ULL;
}

}  // namespace mapper_speed

#endif

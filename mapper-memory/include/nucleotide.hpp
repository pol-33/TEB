#ifndef MAPPER_MEMORY_NUCLEOTIDE_HPP
#define MAPPER_MEMORY_NUCLEOTIDE_HPP

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace mapper_memory {

constexpr uint8_t kSentinelRank = 0;
constexpr uint8_t kSeparatorRank = 1;
constexpr uint8_t kARank = 2;
constexpr uint8_t kCRank = 3;
constexpr uint8_t kGRank = 4;
constexpr uint8_t kTRank = 5;
constexpr uint8_t kAlphabetSize = 6;

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

inline uint8_t rank_from_base(char base) {
    switch (normalize_base(base)) {
        case 'A':
            return kARank;
        case 'C':
            return kCRank;
        case 'G':
            return kGRank;
        case 'T':
            return kTRank;
        default:
            return kARank;
    }
}

inline uint8_t packed_from_rank(uint8_t rank) {
    if (rank < kARank || rank > kTRank) {
        throw std::runtime_error("packed_from_rank received a non-DNA rank");
    }
    return static_cast<uint8_t>(rank - kARank);
}

inline uint8_t rank_from_packed(uint8_t packed) {
    return static_cast<uint8_t>(packed + kARank);
}

inline char base_from_rank(uint8_t rank) {
    switch (rank) {
        case kARank:
            return 'A';
        case kCRank:
            return 'C';
        case kGRank:
            return 'G';
        case kTRank:
            return 'T';
        case kSeparatorRank:
            return '\1';
        case kSentinelRank:
            return '\0';
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

inline std::size_t packed_byte_count(uint64_t symbols) {
    return static_cast<std::size_t>((symbols + 3ULL) / 4ULL);
}

inline void ensure_packed_capacity(std::vector<uint8_t>& packed, uint64_t symbol_count) {
    const std::size_t required = packed_byte_count(symbol_count);
    if (packed.size() < required) {
        packed.resize(required, 0);
    }
}

inline void set_packed_rank(std::vector<uint8_t>& packed, uint64_t index, uint8_t rank) {
    ensure_packed_capacity(packed, index + 1ULL);
    const std::size_t byte_index = static_cast<std::size_t>(index / 4ULL);
    const uint8_t shift = static_cast<uint8_t>((index % 4ULL) * 2ULL);
    const uint8_t value = packed_from_rank(rank);
    packed[byte_index] &= static_cast<uint8_t>(~(0x3u << shift));
    packed[byte_index] |= static_cast<uint8_t>(value << shift);
}

inline uint8_t get_packed_rank(const uint8_t* packed, uint64_t index) {
    const std::size_t byte_index = static_cast<std::size_t>(index / 4ULL);
    const uint8_t shift = static_cast<uint8_t>((index % 4ULL) * 2ULL);
    const uint8_t packed_value = static_cast<uint8_t>((packed[byte_index] >> shift) & 0x3u);
    return rank_from_packed(packed_value);
}

inline void append_packed_rank(std::vector<uint8_t>& packed, uint64_t index, uint8_t rank) {
    set_packed_rank(packed, index, rank);
}

inline std::vector<uint8_t> pack_bases(const std::string& sequence) {
    std::vector<uint8_t> packed;
    packed.reserve(packed_byte_count(sequence.size()));
    for (uint64_t i = 0; i < sequence.size(); ++i) {
        append_packed_rank(packed, i, rank_from_base(sequence[static_cast<std::size_t>(i)]));
    }
    return packed;
}

}  // namespace mapper_memory

#endif

#ifndef MAPPER_MEMORY_NUCLEOTIDE_HPP
#define MAPPER_MEMORY_NUCLEOTIDE_HPP

#include <cctype>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

// ============================================================================
// Bit-manipulation tricks for nucleotide processing
// Based on ASCII encoding properties:
//   A -> 0x41 (0100 0001)  a -> 0x61
//   C -> 0x43 (0100 0011)  c -> 0x63  
//   G -> 0x47 (0100 0111)  g -> 0x67
//   T -> 0x54 (0101 0100)  t -> 0x74
//   N -> 0x4E (0100 1110)  n -> 0x6E
// ============================================================================

// Fast uppercase: clear bit 5 (0x20) to convert lowercase to uppercase
// Works because lowercase letters differ from uppercase only in bit 5
inline constexpr char fast_toupper(char c) {
    return static_cast<char>(c & ~0x20);
}

// Lookup table for base -> rank conversion (pre-computed for speed)
// Index by (uppercase_char & 0x1F) to get a 5-bit index
alignas(64) inline constexpr uint8_t kBaseToRank[32] = {
    kARank, kARank, kARank, kARank, kARank, kARank, kARank, kARank,  // 0-7
    kARank, kARank, kARank, kARank, kARank, kARank, kARank, kARank,  // 8-15
    kARank, kARank, kARank, kARank, kTRank, kARank, kARank, kARank,  // 16-23 (T=20->kTRank)
    kARank, kARank, kARank, kARank, kARank, kARank, kARank, kARank   // 24-31
};

// Direct lookup: A=1, C=3, G=7, T=20 (after &0x1F)
// Build proper lookup indexed by (char & 0x1F)
// A&0x1F=1, C&0x1F=3, G&0x1F=7, T&0x1F=20, N&0x1F=14
alignas(64) inline constexpr uint8_t kCharToRank[32] = {
    kARank, kARank, kARank, kCRank, kARank, kARank, kARank, kGRank,  // 0-7: C=3, G=7
    kARank, kARank, kARank, kARank, kARank, kARank, kARank, kARank,  // 8-15: N=14->A
    kARank, kARank, kARank, kARank, kTRank, kARank, kARank, kARank,  // 16-23: T=20
    kARank, kARank, kARank, kARank, kARank, kARank, kARank, kARank   // 24-31
};

// Fast rank lookup using bit manipulation
// Works for both upper and lowercase
inline uint8_t rank_from_base_fast(char base) {
    // Clear bit 5 for case-insensitivity, then use bottom 5 bits as index
    return kCharToRank[(base & ~0x20) & 0x1F];
}

// Lookup table for fast normalization: char -> normalized uppercase ACGT
// Index by (char & 0x1F), maps to 'A', 'C', 'G', or 'T'
alignas(64) inline constexpr char kNormalizeBase[32] = {
    'A', 'A', 'A', 'C', 'A', 'A', 'A', 'G',  // 0-7: A=1->A, C=3->C, G=7->G
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',  // 8-15: N=14->A
    'A', 'A', 'A', 'A', 'T', 'A', 'A', 'A',  // 16-23: T=20->T
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'   // 24-31
};

// Fast normalize using lookup table
inline char normalize_base_fast(char base) {
    return kNormalizeBase[(base & ~0x20) & 0x1F];
}

// Original normalize_base for compatibility
inline char normalize_base(char base) {
    return normalize_base_fast(base);
}

// Use fast lookup for rank_from_base
inline uint8_t rank_from_base(char base) {
    return rank_from_base_fast(base);
}

// Fast complement lookup table
// A<->T, C<->G (only works on normalized uppercase ACGT)
// A=0x41, C=0x43, G=0x47, T=0x54
// Complement: A->T, C->G, G->C, T->A
alignas(64) inline constexpr char kComplementBase[32] = {
    'A', 'T', 'A', 'G', 'A', 'A', 'A', 'C',  // 0-7: A=1->T, C=3->G, G=7->C
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',  // 8-15
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A',  // 16-23: T=20->A
    'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'   // 24-31
};

// Fast complement for normalized base (uppercase ACGT only)
inline char complement_base_fast(char base) {
    // For normalized input: A(0x41)&0x1F=1, C(0x43)&0x1F=3, G(0x47)&0x1F=7, T(0x54)&0x1F=20
    return kComplementBase[base & 0x1F];
}

// Fast reverse complement for normalized sequence
inline std::string reverse_complement_fast(const std::string& normalized) {
    const std::size_t len = normalized.size();
    std::string out(len, 'A');
    const char* src = normalized.data();
    char* dst = &out[0];
    
    for (std::size_t i = 0; i < len; ++i) {
        dst[i] = kComplementBase[src[len - 1 - i] & 0x1F];
    }
    return out;
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

// Fast bulk sequence normalization
// Uses lookup table instead of switch statements
inline std::string normalize_sequence(const std::string& sequence) {
    std::string out;
    out.resize(sequence.size());
    const char* src = sequence.data();
    char* dst = &out[0];
    const std::size_t len = sequence.size();
    
    // Process in bulk using the lookup table
    for (std::size_t i = 0; i < len; ++i) {
        dst[i] = kNormalizeBase[(src[i] & ~0x20) & 0x1F];
    }
    return out;
}

// In-place normalization for pre-allocated buffer
inline void normalize_sequence_inplace(char* data, std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) {
        data[i] = kNormalizeBase[(data[i] & ~0x20) & 0x1F];
    }
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

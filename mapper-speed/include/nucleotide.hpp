#ifndef MAPPER_SPEED_NUCLEOTIDE_HPP
#define MAPPER_SPEED_NUCLEOTIDE_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

inline const uint8_t* ascii_to_code_table() {
    static uint8_t table[256] = {};
    static bool initialized = false;
    if (!initialized) {
        std::fill(std::begin(table), std::end(table), kBaseCodeInvalid);
        for (int c : {'A', 'a'}) table[static_cast<uint8_t>(c)] = 0;
        for (int c : {'C', 'c'}) table[static_cast<uint8_t>(c)] = 1;
        for (int c : {'G', 'g'}) table[static_cast<uint8_t>(c)] = 2;
        for (int c : {'T', 't'}) table[static_cast<uint8_t>(c)] = 3;
        initialized = true;
    }
    return table;
}

inline const char* normalize_ascii_table() {
    static char table[256] = {};
    static bool initialized = false;
    if (!initialized) {
        std::fill(std::begin(table), std::end(table), 'N');
        table[static_cast<uint8_t>('A')] = table[static_cast<uint8_t>('a')] = 'A';
        table[static_cast<uint8_t>('C')] = table[static_cast<uint8_t>('c')] = 'C';
        table[static_cast<uint8_t>('G')] = table[static_cast<uint8_t>('g')] = 'G';
        table[static_cast<uint8_t>('T')] = table[static_cast<uint8_t>('t')] = 'T';
        initialized = true;
    }
    return table;
}

inline const char* complement_ascii_table() {
    static char table[256] = {};
    static bool initialized = false;
    if (!initialized) {
        std::fill(std::begin(table), std::end(table), 'N');
        table[static_cast<uint8_t>('A')] = table[static_cast<uint8_t>('a')] = 'T';
        table[static_cast<uint8_t>('C')] = table[static_cast<uint8_t>('c')] = 'G';
        table[static_cast<uint8_t>('G')] = table[static_cast<uint8_t>('g')] = 'C';
        table[static_cast<uint8_t>('T')] = table[static_cast<uint8_t>('t')] = 'A';
        initialized = true;
    }
    return table;
}

inline uint8_t base_to_code(char base) {
    return ascii_to_code_table()[static_cast<uint8_t>(base)];
}

inline char normalize_base(char base) {
    return normalize_ascii_table()[static_cast<uint8_t>(base)];
}

inline char complement_base(char base) {
    return complement_ascii_table()[static_cast<uint8_t>(base)];
}

inline std::string normalize_sequence(const std::string& sequence) {
    std::string out(sequence.size(), 'N');
    const char* table = normalize_ascii_table();
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        out[i] = table[static_cast<uint8_t>(sequence[i])];
    }
    return out;
}

inline std::string reverse_complement(const std::string& sequence) {
    std::string out(sequence.size(), 'N');
    const char* table = complement_ascii_table();
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        out[i] = table[static_cast<uint8_t>(sequence[sequence.size() - 1 - i])];
    }
    return out;
}

inline void normalize_and_reverse_complement(const std::string& sequence,
                                             std::string& normalized,
                                             std::string& reverse_comp) {
    normalized.resize(sequence.size());
    reverse_comp.resize(sequence.size());
    const char* normalize_table = normalize_ascii_table();
    const char* complement_table = complement_ascii_table();
    const std::size_t n = sequence.size();
    for (std::size_t i = 0; i < n; ++i) {
        const uint8_t front = static_cast<uint8_t>(sequence[i]);
        const uint8_t back = static_cast<uint8_t>(sequence[n - 1u - i]);
        normalized[i] = normalize_table[front];
        reverse_comp[i] = complement_table[back];
    }
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

inline bool pack_normalized_sequence(const std::string& sequence,
                                     std::vector<uint8_t>& packed) {
    packed.assign(packed_base_bytes(sequence.size()), 0);
    bool all_acgt = true;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const uint8_t code = base_to_code(sequence[i]);
        if (code == kBaseCodeInvalid) {
            all_acgt = false;
            continue;
        }
        packed_set(packed, i, code);
    }
    return all_acgt;
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

#ifndef __UTILS__
#define __UTILS__

#include <unordered_map>
#include <string>
#include <vector>
#include <cstdint>
#include <cstring>

using namespace std;

#define CACHE_BLOCK_SIZE 0
#define NUMBER_CORES 0
#define IO_BUFFER_SIZE (1 << 20)

typedef unordered_map<string,vector<int>> kmer_table_t;

static inline void update_kmer_table(const string& text, kmer_table_t& text_kmer_table, const unsigned int kmer_length) {
    if (text.length() < kmer_length) return;

    string kmer(text, 0, kmer_length);
    text_kmer_table[kmer].push_back(0);

    for (size_t i = 1; i <= text.length() - kmer_length; ++i) {
        kmer.assign(text, i, kmer_length);
        text_kmer_table[kmer].push_back(i);
    }
}

void print_kmer_table(const kmer_table_t& table);

//A -> 65  -> 0b01000001 -> 0x41
//C -> 67  -> 0b01000011 -> 0x43
//G -> 71  -> 0b01000111 -> 0x47
//T -> 84  -> 0b01010100 -> 0x54
//N -> 78  -> 0b01001110 -> 0x4E

//a -> 97  -> 0b01100001 -> 0x61
//c -> 99  -> 0b01100011 -> 0x63
//g -> 103 -> 0b01100111 -> 0x67
//t -> 116 -> 0b01110100 -> 0x74
//n -> 110 -> 0b01101110 -> 0x6E

// bits[1:0]==0b11 uniquely identify G/C (both upper and lower case)
#define gc_matching(char_read_) \
    ((char_read_ & 0x3) == 0x3)

// Bulk GC counter: processes 8 bases per CPU cycle using 64-bit word tricks.
// For each byte, G(0x47) and C(0x43) have bits[1:0]==0b11; A, T, N do not.
// Strategy: extract bit0 and bit1 of every byte independently, AND them,
// then popcount — each byte with both bits set contributes exactly 1.
static inline long long count_gc_bulk(const char* __restrict__ data, size_t len) {
    long long count = 0;
    size_t i = 0;
    // 8-byte vectorised loop
    for (; i + 8 <= len; i += 8) {
        uint64_t chunk;
        memcpy(&chunk, data + i, 8);                               // safe unaligned load
        uint64_t lo =  chunk        & UINT64_C(0x0101010101010101); // bit0 of each byte
        uint64_t hi = (chunk >> 1)  & UINT64_C(0x0101010101010101); // bit1 → bit0 of each byte
        count += __builtin_popcountll(lo & hi);                    // bytes where both bits [1:0] are set
    }
    // Scalar tail for remaining bytes
    for (; i < len; ++i)
        count += gc_matching(data[i]);
    return count;
}

#endif
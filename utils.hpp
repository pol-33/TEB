#ifndef __UTILS__
#define __UTILS__

#include <unordered_map>
#include <string>
#include <list>

using namespace std;

#define CACHE_BLOCK_SIZE 0

#define NUMBER_CORES 0

typedef unordered_map<string,list<int>> kmer_table_t;

static inline void update_kmer_table(const string& text, kmer_table_t& text_kmer_table, const int kmer_length) {
    if (text.length() < kmer_length) return;

    string kmer = text.substr(0, kmer_length);
    text_kmer_table[kmer].push_back(0);

    for (size_t i = 1; i <= text.length() - kmer_length; ++i) {
        // slide window: remove first char, add next char
        kmer.erase(0, 1);
        kmer.push_back(text[i + kmer_length - 1]);

        text_kmer_table[kmer].push_back(i);
    }
}

void print_kmer_table(const kmer_table_t& table);

//A -> 65  -> 0b1000001 -> 0x41
//C -> 67  -> 0b1000011 -> 0x43
//G -> 71  -> 0b1000111 -> 0x47
//T -> 84  -> 0b1010100 -> 0x54
//N -> 78  -> 0b1001110 -> 0x4E

//a -> 97  -> 0b1100001 -> 0x61
//c -> 99  -> 0b1100011 -> 0x63
//g -> 103 -> 0b1100111 -> 0x67
//t -> 116 -> 0b1101110 -> 0x74
//n -> 110 -> 0b1110100 -> 0x6E

#define gc_matching(char_read_) \
    ((char_read_ & 0x3) == 0x3)

#endif
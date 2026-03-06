#ifndef __COMMON__
#define __COMMON__

#include <string>
#include <unordered_map>

struct SequenceStats {
    long long length;
    long long gc_count;
    double gc_content;
};

struct GlobalStats {
    int num_sequences;
    long long total_length;
    long long total_gc_count;
    double overall_gc_content;
    std::unordered_map<std::string, SequenceStats> per_sequence;
};

#endif

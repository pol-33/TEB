#ifndef MAPPER_MEMORY_BACKWARD_SEARCH_HPP
#define MAPPER_MEMORY_BACKWARD_SEARCH_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "fm_index.hpp"

namespace mapper_memory {

struct SearchResult {
    uint64_t sa_lo = 0;
    uint64_t sa_hi = 0;
    int edit_dist = 0;
    std::string cigar_ops;
};

std::pair<uint64_t, uint64_t> exact_search(const FMIndexView& index, const std::string& read);

void inexact_search(const FMIndexView& index,
                    const std::string& read,
                    int max_errors,
                    std::vector<SearchResult>& results);

}  // namespace mapper_memory

#endif

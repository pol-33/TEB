#ifndef MAPPER_MEMORY_FM_SEARCH_HPP
#define MAPPER_MEMORY_FM_SEARCH_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "fm_index.hpp"

namespace mapper_memory {

struct SAInterval {
    uint64_t lo = 0;
    uint64_t hi = 0;
};

struct SearchResult {
    uint64_t sa_lo = 0;
    uint64_t sa_hi = 0;
    int edit_dist = 0;
};

SAInterval exact_search(const FMIndexView::ChromosomeView& index, const std::string& read);
void inexact_search(const FMIndexView::ChromosomeView& index,
                    const std::string& read,
                    int max_errors,
                    std::vector<SearchResult>& results,
                    std::size_t max_results = 64U);

}  // namespace mapper_memory

#endif

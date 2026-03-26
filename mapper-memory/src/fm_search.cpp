#include "fm_search.hpp"

#include <algorithm>
#include <string_view>
#include <vector>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

SAInterval exact_search_range(const FMIndexView::ChromosomeView& index, std::string_view read) {
    uint64_t lo = 0;
    uint64_t hi = index.n();

    for (std::size_t i = read.size(); i > 0 && lo < hi; --i) {
        const uint8_t rank = rank_from_base(read[i - 1U]);
        lo = index.c_array[rank] + index.occ(rank, lo);
        hi = index.c_array[rank] + index.occ(rank, hi);
    }
    return {lo, hi};
}

std::vector<int> compute_lower_bounds(const FMIndexView::ChromosomeView& index, const std::string& read) {
    std::vector<int> lower_bounds(read.size(), 0);
    std::size_t segment_start = 0;
    int segments = 0;

    for (std::size_t i = 0; i < read.size(); ++i) {
        const SAInterval interval = exact_search_range(
            index,
            std::string_view(read.data() + segment_start, i - segment_start + 1U));
        if (interval.lo == interval.hi) {
            ++segments;
            segment_start = i + 1U;
        }
        lower_bounds[i] = segments;
    }
    return lower_bounds;
}

void maybe_add_result(std::vector<SearchResult>& results,
                      uint64_t lo,
                      uint64_t hi,
                      int edit_dist) {
    for (SearchResult& result : results) {
        if (result.sa_lo == lo && result.sa_hi == hi) {
            result.edit_dist = std::min(result.edit_dist, edit_dist);
            return;
        }
    }
    results.push_back(SearchResult{lo, hi, edit_dist});
}

void recurse_search(const FMIndexView::ChromosomeView& index,
                    const std::string& read,
                    int pos,
                    int edits_used,
                    int max_errors,
                    uint64_t lo,
                    uint64_t hi,
                    const std::vector<int>& lower_bounds,
                    std::vector<SearchResult>& results,
                    std::size_t max_results) {
    if (lo >= hi || edits_used > max_errors || results.size() >= max_results) {
        return;
    }
    if (pos >= 0 && edits_used + lower_bounds[static_cast<std::size_t>(pos)] > max_errors) {
        return;
    }
    if (pos < 0) {
        maybe_add_result(results, lo, hi, edits_used);
        return;
    }

    const uint8_t target = rank_from_base(read[static_cast<std::size_t>(pos)]);
    auto step = [&](uint8_t rank, int next_pos, int cost) {
        const uint64_t next_lo = index.c_array[rank] + index.occ(rank, lo);
        const uint64_t next_hi = index.c_array[rank] + index.occ(rank, hi);
        if (next_lo < next_hi) {
            recurse_search(index, read, next_pos, edits_used + cost, max_errors, next_lo, next_hi,
                           lower_bounds, results, max_results);
        }
    };

    step(target, pos - 1, 0);
    if (edits_used == max_errors) {
        return;
    }
    for (uint8_t rank = kARank; rank <= kTRank; ++rank) {
        if (rank != target) {
            step(rank, pos - 1, 1);
        }
    }
    recurse_search(index, read, pos - 1, edits_used + 1, max_errors, lo, hi, lower_bounds, results, max_results);
    for (uint8_t rank = kARank; rank <= kTRank; ++rank) {
        step(rank, pos, 1);
    }
}

}  // namespace

SAInterval exact_search(const FMIndexView::ChromosomeView& index, const std::string& read) {
    return exact_search_range(index, read);
}

void inexact_search(const FMIndexView::ChromosomeView& index,
                    const std::string& read,
                    int max_errors,
                    std::vector<SearchResult>& results,
                    std::size_t max_results) {
    results.clear();
    if (read.empty()) {
        return;
    }
    const std::vector<int> lower_bounds = compute_lower_bounds(index, read);
    recurse_search(index,
                   read,
                   static_cast<int>(read.size()) - 1,
                   0,
                   max_errors,
                   0,
                   index.n(),
                   lower_bounds,
                   results,
                   max_results);

    std::sort(results.begin(), results.end(), [](const SearchResult& lhs, const SearchResult& rhs) {
        if (lhs.edit_dist != rhs.edit_dist) {
            return lhs.edit_dist < rhs.edit_dist;
        }
        if (lhs.sa_lo != rhs.sa_lo) {
            return lhs.sa_lo < rhs.sa_lo;
        }
        return lhs.sa_hi < rhs.sa_hi;
    });
    if (results.size() > max_results) {
        results.resize(max_results);
    }
}

}  // namespace mapper_memory

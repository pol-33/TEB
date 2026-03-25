#include "backward_search.hpp"

#include <algorithm>
#include <array>
#include <string_view>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

std::pair<uint64_t, uint64_t> exact_search_range(const FMIndexView& index, std::string_view read) {
    uint64_t lo = 0;
    uint64_t hi = index.n();
    for (std::size_t i = read.size(); i > 0; --i) {
        const uint8_t code = encode_base(read[i - 1]);
        lo = static_cast<uint64_t>(index.c_array()[code]) + index.occ(code, lo);
        hi = static_cast<uint64_t>(index.c_array()[code]) + index.occ(code, hi);
        if (lo >= hi) {
            return {0, 0};
        }
    }
    return {lo, hi};
}

std::vector<int> compute_lower_bounds(const FMIndexView& index, const std::string& read) {
    std::vector<int> lower_bounds(read.size(), 0);
    int segments = 0;
    std::size_t segment_start = 0;

    for (std::size_t end = 0; end < read.size(); ++end) {
        const auto range = exact_search_range(index, std::string_view(read.data() + segment_start, end - segment_start + 1));
        if (range.first == range.second) {
            ++segments;
            segment_start = end + 1;
        }
        lower_bounds[end] = segments;
    }

    return lower_bounds;
}

void maybe_add_result(std::vector<SearchResult>& results, uint64_t lo, uint64_t hi, int edits, const std::string& ops) {
    for (SearchResult& result : results) {
        if (result.sa_lo == lo && result.sa_hi == hi) {
            if (edits < result.edit_dist) {
                result.edit_dist = edits;
                result.cigar_ops = ops;
            }
            return;
        }
    }
    results.push_back(SearchResult{lo, hi, edits, ops});
}

void recurse_search(const FMIndexView& index,
                    const std::string& read,
                    int pos,
                    uint64_t lo,
                    uint64_t hi,
                    int edits,
                    int max_errors,
                    const std::vector<int>& lower_bounds,
                    std::string& ops,
                    std::vector<SearchResult>& results) {
    if (lo >= hi || edits > max_errors || results.size() >= 32U) {
        return;
    }
    if (pos >= 0 && edits + lower_bounds[static_cast<std::size_t>(pos)] > max_errors) {
        return;
    }
    if (pos < 0) {
        std::string forward_ops(ops.rbegin(), ops.rend());
        maybe_add_result(results, lo, hi, edits, forward_ops);
        return;
    }

    const uint8_t target = encode_base(read[static_cast<std::size_t>(pos)]);
    auto step = [&](uint8_t code, int next_pos, int edit_penalty, char op) {
        const uint64_t next_lo = static_cast<uint64_t>(index.c_array()[code]) + index.occ(code, lo);
        const uint64_t next_hi = static_cast<uint64_t>(index.c_array()[code]) + index.occ(code, hi);
        if (next_lo >= next_hi) {
            return;
        }
        ops.push_back(op);
        recurse_search(index, read, next_pos, next_lo, next_hi, edits + edit_penalty, max_errors, lower_bounds, ops, results);
        ops.pop_back();
    };

    step(target, pos - 1, 0, 'M');
    if (edits == max_errors) {
        return;
    }

    for (uint8_t code = kACode; code <= kTCode; ++code) {
        if (code == target) {
            continue;
        }
        step(code, pos - 1, 1, 'M');
    }

    ops.push_back('I');
    recurse_search(index, read, pos - 1, lo, hi, edits + 1, max_errors, lower_bounds, ops, results);
    ops.pop_back();

    for (uint8_t code = kACode; code <= kTCode; ++code) {
        step(code, pos, 1, 'D');
    }
}

}  // namespace

std::pair<uint64_t, uint64_t> exact_search(const FMIndexView& index, const std::string& read) {
    return exact_search_range(index, read);
}

void inexact_search(const FMIndexView& index,
                    const std::string& read,
                    int max_errors,
                    std::vector<SearchResult>& results) {
    results.clear();
    if (read.empty()) {
        return;
    }

    const std::vector<int> lower_bounds = compute_lower_bounds(index, read);
    std::string ops;
    ops.reserve(read.size() + static_cast<std::size_t>(max_errors));
    recurse_search(index, read, static_cast<int>(read.size()) - 1, 0, index.n(), 0, max_errors, lower_bounds, ops, results);

    std::sort(results.begin(), results.end(), [](const SearchResult& lhs, const SearchResult& rhs) {
        if (lhs.edit_dist != rhs.edit_dist) {
            return lhs.edit_dist < rhs.edit_dist;
        }
        if (lhs.sa_lo != rhs.sa_lo) {
            return lhs.sa_lo < rhs.sa_lo;
        }
        return lhs.sa_hi < rhs.sa_hi;
    });
    if (results.size() > 32U) {
        results.resize(32U);
    }
}

}  // namespace mapper_memory

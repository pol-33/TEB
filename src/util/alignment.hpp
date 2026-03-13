#pragma once

#include <algorithm>
#include <string>
#include <vector>

struct AlignResult {
    int         edit_dist;
    std::string cigar;
};

namespace detail_alignment {

static inline std::string rle_cigar(const std::vector<char>& ops_rev) {
    if (ops_rev.empty()) return "";

    std::string out;
    char cur = ops_rev.back();
    int run = 0;
    for (auto it = ops_rev.rbegin(); it != ops_rev.rend(); ++it) {
        if (*it == cur) {
            ++run;
        } else {
            out += std::to_string(run);
            out.push_back(cur);
            cur = *it;
            run = 1;
        }
    }
    out += std::to_string(run);
    out.push_back(cur);
    return out;
}

} // namespace detail_alignment

// Banded DP alignment of pattern vs text. Returns k+1 if no valid alignment.
static inline AlignResult align_banded(
    const std::string& pattern,
    const std::string& text,
    int                max_errors
) {
    const int m = static_cast<int>(pattern.size());
    const int n = static_cast<int>(text.size());
    const int INF = max_errors + 1;

    if (max_errors < 0) return {INF, ""};
    if (m == 0) return {0, "0M"};
    if (n == 0) return {(m <= max_errors) ? m : INF, (m <= max_errors) ? std::to_string(m) + "I" : ""};

    const int cols = n + 1;
    const size_t total_cells = static_cast<size_t>(m + 1) * static_cast<size_t>(cols);

    thread_local std::vector<int> dp;
    thread_local std::vector<char> bt;
    dp.assign(total_cells, INF);
    bt.assign(total_cells, 0);

    auto idx = [cols](int i, int j) -> size_t {
        return static_cast<size_t>(i) * static_cast<size_t>(cols) + static_cast<size_t>(j);
    };

    dp[idx(0, 0)] = 0;
    for (int i = 1; i <= m && i <= max_errors; ++i) {
        dp[idx(i, 0)] = i;
        bt[idx(i, 0)] = 'I';
    }
    for (int j = 1; j <= n && j <= max_errors; ++j) {
        dp[idx(0, j)] = j;
        bt[idx(0, j)] = 'D';
    }

    for (int i = 1; i <= m; ++i) {
        const int j_lo = std::max(1, i - max_errors);
        const int j_hi = std::min(n, i + max_errors);
        for (int j = j_lo; j <= j_hi; ++j) {
            int best = INF;
            char op = 0;

            const int sub = dp[idx(i - 1, j - 1)] + ((pattern[i - 1] == text[j - 1]) ? 0 : 1);
            if (sub < best) {
                best = sub;
                op = 'M';
            }

            const int ins = dp[idx(i - 1, j)] + 1;
            if (ins < best) {
                best = ins;
                op = 'I';
            }

            const int del = dp[idx(i, j - 1)] + 1;
            if (del < best) {
                best = del;
                op = 'D';
            }

            dp[idx(i, j)] = best;
            bt[idx(i, j)] = op;
        }
    }

    int best_j = -1;
    int best_ed = INF;
    for (int j = std::max(0, m - max_errors); j <= std::min(n, m + max_errors); ++j) {
        if (dp[idx(m, j)] < best_ed) {
            best_ed = dp[idx(m, j)];
            best_j = j;
        }
    }
    if (best_j < 0 || best_ed > max_errors) return {INF, ""};

    std::vector<char> ops_rev;
    int i = m;
    int j = best_j;
    while (i > 0 || j > 0) {
        char op = bt[idx(i, j)];
        if (op == 0) break;
        ops_rev.push_back(op);
        if (op == 'M') {
            --i; --j;
        } else if (op == 'I') {
            --i;
        } else {
            --j;
        }
    }

    return {best_ed, detail_alignment::rle_cigar(ops_rev)};
}

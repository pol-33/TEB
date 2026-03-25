#include "alignment.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace mapper_memory {

namespace {

constexpr int kInf = std::numeric_limits<int>::max() / 4;

std::string compress_cigar(const std::string& ops) {
    if (ops.empty()) {
        return "*";
    }
    std::string cigar;
    int count = 1;
    for (std::size_t i = 1; i <= ops.size(); ++i) {
        if (i < ops.size() && ops[i] == ops[i - 1]) {
            ++count;
            continue;
        }
        cigar += std::to_string(count);
        cigar.push_back(ops[i - 1]);
        count = 1;
    }
    return cigar;
}

}  // namespace

void AlignmentWorkspace::reserve(std::size_t read_length, std::size_t ref_length) {
    if (read_length + 1 <= rows_ && ref_length + 1 <= cols_) {
        return;
    }
    rows_ = std::max(rows_, read_length + 1);
    cols_ = std::max(cols_, ref_length + 1);
    dp_.assign(rows_ * cols_, kInf);
    trace_.assign(rows_ * cols_, 0);
}

Alignment AlignmentWorkspace::align(const std::string& read, const std::string& ref_window, int max_errors) {
    reserve(read.size(), ref_window.size());
    const std::size_t n = read.size();
    const std::size_t m = ref_window.size();
    const std::size_t matrix_size = rows_ * cols_;
    std::fill(dp_.begin(), dp_.begin() + matrix_size, kInf);
    std::fill(trace_.begin(), trace_.begin() + matrix_size, 0);

    auto cell = [&](std::size_t i, std::size_t j) -> int& {
        return dp_[i * cols_ + j];
    };
    auto trace = [&](std::size_t i, std::size_t j) -> uint8_t& {
        return trace_[i * cols_ + j];
    };

    cell(0, 0) = 0;
    for (std::size_t i = 1; i <= n && static_cast<int>(i) <= max_errors; ++i) {
        cell(i, 0) = static_cast<int>(i);
        trace(i, 0) = 'I';
    }
    for (std::size_t j = 1; j <= m && static_cast<int>(j) <= max_errors; ++j) {
        cell(0, j) = static_cast<int>(j);
        trace(0, j) = 'D';
    }

    for (std::size_t i = 1; i <= n; ++i) {
        const std::size_t j_begin = (i > static_cast<std::size_t>(max_errors)) ? i - static_cast<std::size_t>(max_errors) : 1U;
        const std::size_t j_end = std::min(m, i + static_cast<std::size_t>(max_errors));
        for (std::size_t j = j_begin; j <= j_end; ++j) {
            int best = cell(i - 1, j - 1) + (read[i - 1] == ref_window[j - 1] ? 0 : 1);
            uint8_t best_op = 'M';

            if (cell(i - 1, j) + 1 < best) {
                best = cell(i - 1, j) + 1;
                best_op = 'I';
            }
            if (cell(i, j - 1) + 1 < best) {
                best = cell(i, j - 1) + 1;
                best_op = 'D';
            }

            cell(i, j) = best;
            trace(i, j) = best_op;
        }
    }

    Alignment result;
    result.edit_dist = cell(n, m);
    if (result.edit_dist >= kInf) {
        throw std::runtime_error("banded alignment failed to find a valid path");
    }

    std::string ops;
    ops.reserve(n + m);
    std::size_t i = n;
    std::size_t j = m;
    while (i > 0 || j > 0) {
        const uint8_t op = trace(i, j);
        if (op == 'M') {
            ops.push_back('M');
            --i;
            --j;
        } else if (op == 'I') {
            ops.push_back('I');
            --i;
        } else if (op == 'D') {
            ops.push_back('D');
            --j;
        } else {
            throw std::runtime_error("invalid backtrace in banded alignment");
        }
    }
    std::reverse(ops.begin(), ops.end());
    result.cigar = compress_cigar(ops);
    return result;
}

Alignment band_align(const std::string& read,
                     const std::string& ref_window,
                     int max_errors,
                     AlignmentWorkspace& workspace) {
    return workspace.align(read, ref_window, max_errors);
}

}  // namespace mapper_memory

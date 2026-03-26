#include "bwt.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mapper_memory {

namespace {

constexpr uint32_t kInvalid = std::numeric_limits<uint32_t>::max();

bool is_lms_pos(const std::vector<uint64_t>& bits, uint32_t pos) {
    if (pos == 0) {
        return false;
    }
    const std::size_t word = static_cast<std::size_t>(pos / 64U);
    const uint32_t bit = pos % 64U;
    return ((bits[word] >> bit) & 1ULL) != 0ULL;
}

std::vector<uint32_t> build_lms_superblocks(const std::vector<uint64_t>& bits) {
    std::vector<uint32_t> superblocks((bits.size() + 7U) / 8U + 1U, 0U);
    uint32_t running = 0;
    for (std::size_t i = 0; i < bits.size(); ++i) {
        if (i % 8U == 0U) {
            superblocks[i / 8U] = running;
        }
        running += static_cast<uint32_t>(__builtin_popcountll(bits[i]));
    }
    superblocks.back() = running;
    return superblocks;
}

uint32_t rank_lms(const std::vector<uint64_t>& bits,
                  const std::vector<uint32_t>& superblocks,
                  uint32_t pos) {
    const std::size_t word = static_cast<std::size_t>(pos / 64U);
    const uint32_t bit = pos % 64U;
    const std::size_t super = word / 8U;
    uint32_t count = superblocks[super];
    for (std::size_t w = super * 8U; w < word; ++w) {
        count += static_cast<uint32_t>(__builtin_popcountll(bits[w]));
    }
    if (bit != 0U) {
        count += static_cast<uint32_t>(__builtin_popcountll(bits[word] & ((1ULL << bit) - 1ULL)));
    }
    return count;
}

std::vector<uint32_t> sa_naive(const std::vector<uint32_t>& s) {
    const uint32_t n = static_cast<uint32_t>(s.size());
    std::vector<uint32_t> sa(n);
    for (uint32_t i = 0; i < n; ++i) {
        sa[i] = i;
    }

    std::sort(sa.begin(), sa.end(), [&](uint32_t lhs, uint32_t rhs) {
        if (lhs == rhs) {
            return false;
        }
        while (lhs < n && rhs < n) {
            if (s[lhs] != s[rhs]) {
                return s[lhs] < s[rhs];
            }
            ++lhs;
            ++rhs;
        }
        return lhs == n;
    });
    return sa;
}

std::vector<uint32_t> sa_doubling(const std::vector<uint32_t>& s) {
    const uint32_t n = static_cast<uint32_t>(s.size());
    std::vector<uint32_t> sa(n);
    std::vector<int64_t> rank(n);
    std::vector<int64_t> tmp(n, 0);

    for (uint32_t i = 0; i < n; ++i) {
        sa[i] = i;
        rank[i] = static_cast<int64_t>(s[i]);
    }

    for (uint32_t k = 1;; k <<= 1U) {
        auto cmp = [&](uint32_t lhs, uint32_t rhs) {
            if (rank[lhs] != rank[rhs]) {
                return rank[lhs] < rank[rhs];
            }
            const int64_t lhs_next = (lhs + k < n) ? rank[lhs + k] : -1;
            const int64_t rhs_next = (rhs + k < n) ? rank[rhs + k] : -1;
            return lhs_next < rhs_next;
        };
        std::sort(sa.begin(), sa.end(), cmp);

        tmp[sa[0]] = 0;
        for (uint32_t i = 1; i < n; ++i) {
            tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        }
        rank.swap(tmp);
        if (rank[sa[n - 1]] == static_cast<int64_t>(n - 1U)) {
            break;
        }
    }

    return sa;
}

std::vector<uint32_t> sa_is_impl(const std::vector<uint32_t>& s, uint32_t upper) {
    const uint32_t n = static_cast<uint32_t>(s.size());
    if (n == 0U) {
        return {};
    }
    if (n == 1U) {
        return {0U};
    }
    if (n == 2U) {
        return (s[0] < s[1]) ? std::vector<uint32_t>{0U, 1U} : std::vector<uint32_t>{1U, 0U};
    }
    if (n < 10U) {
        return sa_naive(s);
    }
    if (n < 40U) {
        return sa_doubling(s);
    }

    std::vector<uint32_t> sa(n, kInvalid);
    std::vector<uint8_t> ls(n, 0);
    ls[n - 1U] = 1U;
    for (uint64_t i = static_cast<uint64_t>(n) - 1ULL; i-- > 0ULL;) {
        const uint32_t idx = static_cast<uint32_t>(i);
        ls[idx] = (s[idx] == s[idx + 1U]) ? ls[idx + 1U] : static_cast<uint8_t>(s[idx] < s[idx + 1U]);
    }

    std::vector<uint64_t> sum_l(static_cast<std::size_t>(upper) + 2U, 0ULL);
    std::vector<uint64_t> sum_s(static_cast<std::size_t>(upper) + 2U, 0ULL);
    for (uint32_t i = 0; i < n; ++i) {
        if (ls[i] == 0U) {
            ++sum_s[s[i]];
        } else {
            ++sum_l[s[i] + 1U];
        }
    }
    for (uint32_t i = 0; i <= upper + 1U; ++i) {
        sum_s[i] += sum_l[i];
        if (i != upper + 1U) {
            sum_l[i + 1U] += sum_s[i];
        }
    }

    std::vector<uint64_t> lms_bits((static_cast<uint64_t>(n) + 63ULL) / 64ULL, 0ULL);
    std::vector<uint32_t> lms;
    for (uint32_t i = 1; i < n; ++i) {
        if (ls[i - 1U] == 0U && ls[i] != 0U) {
            lms.push_back(i);
            lms_bits[static_cast<std::size_t>(i / 64U)] |= (1ULL << (i % 64U));
        }
    }
    const std::vector<uint32_t> lms_superblocks = build_lms_superblocks(lms_bits);

    auto induce = [&](const std::vector<uint32_t>& ordered_lms) {
        std::fill(sa.begin(), sa.end(), kInvalid);
        std::vector<uint64_t> buf = sum_s;

        for (uint32_t pos : ordered_lms) {
            const uint64_t slot = buf[s[pos]]++;
            sa[static_cast<std::size_t>(slot)] = pos;
        }

        buf = sum_l;
        sa[static_cast<std::size_t>(buf[s[n - 1U]]++)] = n - 1U;
        for (uint32_t i = 0; i < n; ++i) {
            const uint32_t v = sa[i];
            if (v != kInvalid && v > 0U && ls[v - 1U] == 0U) {
                sa[static_cast<std::size_t>(buf[s[v - 1U]]++)] = v - 1U;
            }
        }

        buf = sum_l;
        for (uint64_t i = static_cast<uint64_t>(n); i-- > 0ULL;) {
            const uint32_t v = sa[static_cast<std::size_t>(i)];
            if (v != kInvalid && v > 0U && ls[v - 1U] != 0U) {
                const uint64_t slot = --buf[s[v - 1U] + 1U];
                sa[static_cast<std::size_t>(slot)] = v - 1U;
            }
        }
    };

    induce(lms);

    if (!lms.empty()) {
        std::vector<uint32_t> sorted_lms;
        sorted_lms.reserve(lms.size());
        for (uint32_t v : sa) {
            if (v != kInvalid && is_lms_pos(lms_bits, v)) {
                sorted_lms.push_back(v);
            }
        }

        std::vector<uint32_t> rec_s(lms.size(), 0U);
        uint32_t rec_upper = 0U;
        rec_s[rank_lms(lms_bits, lms_superblocks, sorted_lms[0])] = 0U;

        for (std::size_t i = 1; i < sorted_lms.size(); ++i) {
            const uint32_t lhs = sorted_lms[i - 1U];
            const uint32_t rhs = sorted_lms[i];
            const uint32_t lhs_idx = rank_lms(lms_bits, lms_superblocks, lhs);
            const uint32_t rhs_idx = rank_lms(lms_bits, lms_superblocks, rhs);
            const uint32_t lhs_end = (lhs_idx + 1U < lms.size()) ? lms[lhs_idx + 1U] : n;
            const uint32_t rhs_end = (rhs_idx + 1U < lms.size()) ? lms[rhs_idx + 1U] : n;

            bool same = (lhs_end - lhs) == (rhs_end - rhs);
            uint32_t l = lhs;
            uint32_t r = rhs;
            while (same && l < lhs_end) {
                if (s[l] != s[r] || ls[l] != ls[r]) {
                    same = false;
                    break;
                }
                ++l;
                ++r;
            }
            if (!same) {
                ++rec_upper;
            }
            rec_s[rhs_idx] = rec_upper;
        }

        const std::vector<uint32_t> rec_sa = sa_is_impl(rec_s, rec_upper);
        std::vector<uint32_t> ordered_lms(lms.size());
        for (std::size_t i = 0; i < rec_sa.size(); ++i) {
            ordered_lms[i] = lms[rec_sa[i]];
        }
        induce(ordered_lms);
    }

    return sa;
}

}  // namespace

std::vector<uint32_t> build_suffix_array_sais(const std::vector<uint8_t>& text) {
    if (text.empty()) {
        return {};
    }
    std::vector<uint32_t> symbols(text.size(), 0U);
    for (std::size_t i = 0; i < text.size(); ++i) {
        symbols[i] = text[i];
    }
    return sa_is_impl(symbols, static_cast<uint32_t>(kAlphabetSize - 1U));
}

BWTData build_bwt(const std::vector<uint8_t>& text, const std::vector<uint32_t>& suffix_array) {
    const uint64_t n = static_cast<uint64_t>(text.size());
    if (static_cast<uint64_t>(suffix_array.size()) != n) {
        throw std::runtime_error("suffix array length does not match text length");
    }

    BWTData data;
    data.text_length = n;
    data.packed_bwt.assign(packed_byte_count(n), 0);

    for (uint64_t row = 0; row < n; ++row) {
        const uint32_t sa_value = suffix_array[static_cast<std::size_t>(row)];
        const uint8_t prev_rank = (sa_value == 0U)
            ? text.back()
            : text[static_cast<std::size_t>(sa_value - 1U)];
        ++data.counts[prev_rank];

        if (prev_rank == kSentinelRank) {
            data.primary_index = row;
            set_packed_rank(data.packed_bwt, row, kARank);
        } else if (prev_rank == kSeparatorRank) {
            data.separator_rows.push_back(static_cast<uint32_t>(row));
            set_packed_rank(data.packed_bwt, row, kARank);
        } else {
            set_packed_rank(data.packed_bwt, row, prev_rank);
        }
    }

    return data;
}

}  // namespace mapper_memory

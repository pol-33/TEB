#include "bwt.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

std::vector<int> sa_naive(const std::vector<int>& s) {
    const int n = static_cast<int>(s.size());
    std::vector<int> sa(n);
    for (int i = 0; i < n; ++i) {
        sa[i] = i;
    }
    std::sort(sa.begin(), sa.end(), [&](int lhs, int rhs) {
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

std::vector<int> sa_doubling(const std::vector<int>& s) {
    const int n = static_cast<int>(s.size());
    std::vector<int> sa(n);
    std::vector<int> rank = s;
    std::vector<int> tmp(n, 0);

    for (int i = 0; i < n; ++i) {
        sa[i] = i;
    }

    for (int k = 1;; k <<= 1) {
        auto cmp = [&](int lhs, int rhs) {
            if (rank[lhs] != rank[rhs]) {
                return rank[lhs] < rank[rhs];
            }
            const int lhs_next = lhs + k < n ? rank[lhs + k] : -1;
            const int rhs_next = rhs + k < n ? rank[rhs + k] : -1;
            return lhs_next < rhs_next;
        };
        std::sort(sa.begin(), sa.end(), cmp);

        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        }
        rank.swap(tmp);
        if (rank[sa[n - 1]] == n - 1) {
            break;
        }
    }

    return sa;
}

std::vector<int> sa_is_impl(const std::vector<int>& s, int upper) {
    const int n = static_cast<int>(s.size());
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return {0};
    }
    if (n == 2) {
        if (s[0] < s[1]) {
            return {0, 1};
        }
        return {1, 0};
    }
    if (n < 10) {
        return sa_naive(s);
    }
    if (n < 40) {
        return sa_doubling(s);
    }

    std::vector<int> sa(n, -1);
    std::vector<bool> ls(n, false);
    ls[n - 1] = true;
    for (int i = n - 2; i >= 0; --i) {
        if (s[i] == s[i + 1]) {
            ls[i] = ls[i + 1];
        } else {
            ls[i] = s[i] < s[i + 1];
        }
    }

    std::vector<int> sum_l(upper + 2, 0);
    std::vector<int> sum_s(upper + 2, 0);
    for (int i = 0; i < n; ++i) {
        if (!ls[i]) {
            ++sum_s[s[i]];
        } else {
            ++sum_l[s[i] + 1];
        }
    }
    for (int i = 0; i <= upper + 1; ++i) {
        sum_s[i] += sum_l[i];
        if (i != upper + 1) {
            sum_l[i + 1] += sum_s[i];
        }
    }

    auto induce = [&](const std::vector<int>& lms) {
        std::fill(sa.begin(), sa.end(), -1);
        std::vector<int> buf(upper + 2, 0);

        std::copy(sum_s.begin(), sum_s.end(), buf.begin());
        for (int pos : lms) {
            if (pos == n) {
                continue;
            }
            sa[buf[s[pos]]++] = pos;
        }

        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        sa[buf[s[n - 1]]++] = n - 1;
        for (int i = 0; i < n; ++i) {
            const int v = sa[i];
            if (v >= 1 && !ls[v - 1]) {
                sa[buf[s[v - 1]]++] = v - 1;
            }
        }

        std::copy(sum_l.begin(), sum_l.end(), buf.begin());
        for (int i = n - 1; i >= 0; --i) {
            const int v = sa[i];
            if (v >= 1 && ls[v - 1]) {
                sa[--buf[s[v - 1] + 1]] = v - 1;
            }
        }
    };

    std::vector<int> lms_map(n + 1, -1);
    int m = 0;
    for (int i = 1; i < n; ++i) {
        if (!ls[i - 1] && ls[i]) {
            lms_map[i] = m++;
        }
    }

    std::vector<int> lms;
    lms.reserve(m);
    for (int i = 1; i < n; ++i) {
        if (!ls[i - 1] && ls[i]) {
            lms.push_back(i);
        }
    }

    induce(lms);

    if (m > 0) {
        std::vector<int> sorted_lms;
        sorted_lms.reserve(m);
        for (int v : sa) {
            if (v >= 0 && lms_map[v] != -1) {
                sorted_lms.push_back(v);
            }
        }

        std::vector<int> rec_s(m, 0);
        int rec_upper = 0;
        rec_s[lms_map[sorted_lms[0]]] = 0;

        for (int i = 1; i < m; ++i) {
            int lhs = sorted_lms[i - 1];
            int rhs = sorted_lms[i];
            const int lhs_end = (lms_map[lhs] + 1 < m) ? lms[lms_map[lhs] + 1] : n;
            const int rhs_end = (lms_map[rhs] + 1 < m) ? lms[lms_map[rhs] + 1] : n;

            bool same = (lhs_end - lhs) == (rhs_end - rhs);
            while (same && lhs < lhs_end) {
                if (s[lhs] != s[rhs] || ls[lhs] != ls[rhs]) {
                    same = false;
                    break;
                }
                ++lhs;
                ++rhs;
            }
            if (!same) {
                ++rec_upper;
            }
            rec_s[lms_map[sorted_lms[i]]] = rec_upper;
        }

        std::vector<int> rec_sa = sa_is_impl(rec_s, rec_upper);
        for (int i = 0; i < m; ++i) {
            sorted_lms[i] = lms[rec_sa[i]];
        }
        induce(sorted_lms);
    }

    return sa;
}

}  // namespace

std::vector<uint32_t> build_suffix_array_sais(const std::string& genome) {
    std::vector<int> text;
    text.reserve(genome.size() + 1);
    for (char base : genome) {
        text.push_back(static_cast<int>(encode_base(base)));
    }
    text.push_back(static_cast<int>(kSentinelCode));

    const std::vector<int> sa = sa_is_impl(text, 4);
    std::vector<uint32_t> out;
    out.reserve(sa.size());
    for (int value : sa) {
        out.push_back(static_cast<uint32_t>(value));
    }
    return out;
}

BWTData build_bwt(const std::string& genome, const std::vector<uint32_t>& suffix_array) {
    const uint64_t text_length = static_cast<uint64_t>(genome.size()) + 1ULL;
    if (suffix_array.size() != text_length) {
        throw std::runtime_error("suffix array length does not match genome length + sentinel");
    }

    BWTData data;
    data.text_length = text_length;
    data.suffix_array = suffix_array;
    data.packed_bwt.assign(packed_byte_count(text_length), 0);

    for (uint64_t i = 0; i < text_length; ++i) {
        const uint32_t sa_value = suffix_array[static_cast<std::size_t>(i)];
        if (sa_value == 0) {
            data.primary_index = i;
            set_packed_code(data.packed_bwt, i, kACode);
            continue;
        }
        const uint64_t prev_pos = static_cast<uint64_t>(sa_value - 1);
        set_packed_code(data.packed_bwt, i, encode_base(genome[static_cast<std::size_t>(prev_pos)]));
    }

    return data;
}

}  // namespace mapper_memory

#include "verifier.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "common.hpp"

namespace mapper_speed {

namespace {

inline constexpr uint8_t kBaseToMaskIndex[32] = {
    0xFF, 0,    0xFF, 1,    0xFF, 0xFF, 0xFF, 2,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 3,    0xFF, 0xFF, 0xFF,
    0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
};

inline uint8_t mask_index_from_base(char base) {
    return kBaseToMaskIndex[(static_cast<uint8_t>(base) & static_cast<uint8_t>(~0x20)) & 0x1F];
}

template <bool UseBuiltinCpuHints>
int myers128_impl(const MyersQuery& query, std::string_view ref, int max_errors) {
    if (query.read_length == 0u) {
        return static_cast<int>(ref.size());
    }
    if (query.read_length > 128u) {
        return max_errors + 1;
    }

    uint64_t pv0 = ~uint64_t{0};
    uint64_t mv0 = 0;
    uint64_t pv1 = query.has_hi ? ~uint64_t{0} : 0;
    uint64_t mv1 = 0;
    int score = static_cast<int>(query.read_length);

    for (std::size_t pos = 0; pos < ref.size(); ++pos) {
        const char base = ref[pos];
        const uint8_t idx = mask_index_from_base(base);
        const uint64_t eq0 = (idx < 4u) ? query.masks_lo[idx] : 0;
        const uint64_t xv0 = eq0 | mv0;
        const uint64_t xh0 = (((eq0 & pv0) + pv0) ^ pv0) | eq0;
        uint64_t ph0 = mv0 | ~(xh0 | pv0);
        uint64_t mh0 = pv0 & xh0;

        const uint64_t carry_ph = (ph0 >> 63u) & 1u;
        const uint64_t carry_mh = (mh0 >> 63u) & 1u;
        ph0 = (ph0 << 1u) | 1u;
        mh0 <<= 1u;
        pv0 = mh0 | ~(xv0 | ph0);
        mv0 = ph0 & xv0;

        if (query.has_hi) {
            const uint64_t eq1 = (idx < 4u) ? query.masks_hi[idx] : 0;
            const uint64_t xv1 = eq1 | mv1;
            const uint64_t xh1 = (((eq1 & pv1) + pv1 + carry_ph) ^ pv1) | eq1;
            uint64_t ph1 = mv1 | ~(xh1 | pv1);
            uint64_t mh1 = pv1 & xh1;
            ph1 = (ph1 << 1u) | carry_ph;
            mh1 = (mh1 << 1u) | carry_mh;
            pv1 = mh1 | ~(xv1 | ph1);
            mv1 = ph1 & xv1;
            if ((ph1 & query.last_bit1) != 0u) {
                ++score;
            } else if ((mh1 & query.last_bit1) != 0u) {
                --score;
            }
        } else {
            if ((ph0 & query.last_bit0) != 0u) {
                ++score;
            } else if ((mh0 & query.last_bit0) != 0u) {
                --score;
            }
        }

        if (UseBuiltinCpuHints) {
            const std::size_t remaining = ref.size() - pos - 1u;
            const int lower_bound = score - static_cast<int>(remaining);
            if (lower_bound > max_errors) {
                return lower_bound;
            }
        }
    }

    return score;
}

int myers_generic(const MyersQuery& query, std::string_view ref, int max_errors) {
    return myers128_impl<false>(query, ref, max_errors);
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("popcnt,bmi2")))
int myers_popcnt(const MyersQuery& query, std::string_view ref, int max_errors) {
    return myers128_impl<true>(query, ref, max_errors);
}

__attribute__((target("avx512f,avx512bw,avx512vl,popcnt,bmi2")))
int myers_avx512(const MyersQuery& query, std::string_view ref, int max_errors) {
    return myers128_impl<true>(query, ref, max_errors);
}
#endif

}  // namespace

MyersDispatch resolve_myers_dispatch() {
    MyersDispatch dispatch{};
    dispatch.selected = &myers_generic;
    dispatch.name = "generic";
#if defined(__x86_64__) || defined(__i386__)
    if (detect_simd_level() == SimdLevel::kAvx512 && __builtin_cpu_supports("bmi2")) {
        dispatch.selected = &myers_avx512;
        dispatch.name = "avx512";
    } else if (detect_simd_level() == SimdLevel::kAvx2 && __builtin_cpu_supports("bmi2")) {
        dispatch.selected = &myers_popcnt;
        dispatch.name = "popcnt+bmi2";
    }
#endif
    return dispatch;
}

MyersQuery build_myers_query(std::string_view read) {
    MyersQuery query{};
    query.read_length = read.size();
    query.has_hi = read.size() > 64u;
    if (read.empty()) {
        return query;
    }

    for (std::size_t i = 0; i < read.size(); ++i) {
        const uint8_t idx = mask_index_from_base(read[i]);
        if (idx >= 4u) {
            continue;
        }
        const uint64_t bit = uint64_t{1} << (i & 63u);
        uint64_t* dst = (i < 64u) ? query.masks_lo : query.masks_hi;
        dst[idx] |= bit;
    }

    query.last_bit0 = uint64_t{1} << ((std::min<std::size_t>(read.size(), 64u) - 1u) & 63u);
    query.last_bit1 = query.has_hi ? (uint64_t{1} << ((read.size() - 65u) & 63u)) : 0;
    return query;
}

int bounded_edit_distance(const MyersDispatch& dispatch,
                          const MyersQuery& query,
                          std::string_view ref,
                          int max_errors) {
    return dispatch.selected(query, ref, max_errors);
}

int banded_score_only(std::string_view read,
                      std::string_view ref,
                      int max_errors,
                      std::vector<int>& prev,
                      std::vector<int>& curr) {
    const std::size_t n = read.size();
    const std::size_t m = ref.size();
    constexpr int kInf = std::numeric_limits<int>::max() / 4;

    prev.assign(m + 1u, kInf);
    curr.assign(m + 1u, kInf);
    prev[0] = 0;
    for (std::size_t j = 1; j <= m && static_cast<int>(j) <= max_errors; ++j) {
        prev[j] = static_cast<int>(j);
    }

    for (std::size_t i = 1; i <= n; ++i) {
        std::fill(curr.begin(), curr.end(), kInf);
        if (static_cast<int>(i) <= max_errors) {
            curr[0] = static_cast<int>(i);
        }
        const std::size_t j_begin =
            (i > static_cast<std::size_t>(max_errors)) ? i - static_cast<std::size_t>(max_errors) : 1u;
        const std::size_t j_end = std::min(m, i + static_cast<std::size_t>(max_errors));
        int row_best = kInf;
        for (std::size_t j = j_begin; j <= j_end; ++j) {
            int best = prev[j - 1u] + (read[i - 1u] == ref[j - 1u] ? 0 : 1);
            best = std::min(best, prev[j] + 1);
            best = std::min(best, curr[j - 1u] + 1);
            curr[j] = best;
            row_best = std::min(row_best, best);
        }
        if (row_best > max_errors) {
            return row_best;
        }
        prev.swap(curr);
    }

    return prev[m];
}

int banded_score_only(std::string_view read,
                      std::string_view ref,
                      int max_errors) {
    std::vector<int> prev;
    std::vector<int> curr;
    return banded_score_only(read, ref, max_errors, prev, curr);
}

}  // namespace mapper_speed

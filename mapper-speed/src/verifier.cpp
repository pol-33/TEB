#include "verifier.hpp"

#include <array>
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

int banded_score_only_scalar_impl(std::string_view read,
                                  std::string_view ref,
                                  int max_errors,
                                  std::vector<int>& prev,
                                  std::vector<int>& curr) {
    const std::size_t n = read.size();
    const std::size_t m = ref.size();
    constexpr int kInf = std::numeric_limits<int>::max() / 4;

    if (prev.size() < m + 1u) {
        prev.resize(m + 1u, kInf);
    }
    if (curr.size() < m + 1u) {
        curr.resize(m + 1u, kInf);
    }
    std::fill(prev.begin(), prev.begin() + static_cast<std::ptrdiff_t>(m + 1u), kInf);
    prev[0] = 0;
    for (std::size_t j = 1; j <= m && static_cast<int>(j) <= max_errors; ++j) {
        prev[j] = static_cast<int>(j);
    }
    std::size_t prev_begin = 0;
    std::size_t prev_end = std::min<std::size_t>(m, static_cast<std::size_t>(max_errors));

    for (std::size_t i = 1; i <= n; ++i) {
        const std::size_t j_begin =
            (i > static_cast<std::size_t>(max_errors)) ? i - static_cast<std::size_t>(max_errors) : 1u;
        const std::size_t j_end = std::min(m, i + static_cast<std::size_t>(max_errors));
        const bool row_has_zero = static_cast<int>(i) <= max_errors;
        const std::size_t curr_begin = row_has_zero ? 0u : j_begin;

        if (row_has_zero) {
            curr[0] = static_cast<int>(i);
        } else {
            curr[j_begin - 1u] = kInf;
        }
        int row_best = kInf;
        for (std::size_t j = j_begin; j <= j_end; ++j) {
            const int diag =
                (j - 1u >= prev_begin && j - 1u <= prev_end) ? prev[j - 1u] : kInf;
            const int up =
                (j >= prev_begin && j <= prev_end) ? prev[j] : kInf;
            int best = diag + (read[i - 1u] == ref[j - 1u] ? 0 : 1);
            best = std::min(best, up + 1);
            best = std::min(best, curr[j - 1u] + 1);
            curr[j] = best;
            row_best = std::min(row_best, best);
        }
        if (row_best > max_errors) {
            return row_best;
        }
        prev.swap(curr);
        prev_begin = curr_begin;
        prev_end = j_end;
    }

    return prev[m];
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("popcnt,bmi2")))
int myers_popcnt(const MyersQuery& query, std::string_view ref, int max_errors) {
    return myers128_impl<true>(query, ref, max_errors);
}

__attribute__((target("avx512f,avx512bw,avx512vl")))
int banded_score_only_avx512_impl(std::string_view read,
                                  std::string_view ref,
                                  int max_errors,
                                  std::vector<int>& prev,
                                  std::vector<int>& curr) {
    const std::size_t n = read.size();
    const std::size_t m = ref.size();
    constexpr int kInf = std::numeric_limits<int>::max() / 4;
    constexpr std::size_t kMaxVectorLanes = 16;

    const std::size_t band_lanes = static_cast<std::size_t>(max_errors * 2 + 1);
    // AVX512 only pays off when the active DP band is wide enough to amortize
    // vector setup costs. For mapper-speed's common k<=3 case, the band width
    // is at most 7, so the scalar/PDEP path remains faster on MN5.
    if (band_lanes > kMaxVectorLanes || band_lanes < 12u || n < 96u || m < 96u) {
        return banded_score_only_scalar_impl(read, ref, max_errors, prev, curr);
    }
    if (prev.size() < m + 1u) {
        prev.resize(m + 1u, kInf);
    }
    if (curr.size() < m + 1u) {
        curr.resize(m + 1u, kInf);
    }
    std::fill(prev.begin(), prev.begin() + static_cast<std::ptrdiff_t>(m + 1u), kInf);
    prev[0] = 0;
    for (std::size_t j = 1; j <= m && static_cast<int>(j) <= max_errors; ++j) {
        prev[j] = static_cast<int>(j);
    }
    std::size_t prev_begin = 0;
    std::size_t prev_end = std::min<std::size_t>(m, static_cast<std::size_t>(max_errors));

    alignas(64) std::array<int, kMaxVectorLanes> diag_buf{};
    alignas(64) std::array<int, kMaxVectorLanes> up_buf{};
    alignas(64) std::array<int, kMaxVectorLanes> base_buf{};
    alignas(64) std::array<char, 64> ref_buf{};

    for (std::size_t i = 1; i <= n; ++i) {
        const std::size_t j_begin =
            (i > static_cast<std::size_t>(max_errors)) ? i - static_cast<std::size_t>(max_errors) : 1u;
        const std::size_t j_end = std::min(m, i + static_cast<std::size_t>(max_errors));
        const bool row_has_zero = static_cast<int>(i) <= max_errors;
        const std::size_t curr_begin = row_has_zero ? 0u : j_begin;
        const std::size_t lane_count = j_end - j_begin + 1u;

        if (row_has_zero) {
            curr[0] = static_cast<int>(i);
        } else {
            curr[j_begin - 1u] = kInf;
        }

        for (std::size_t lane = 0; lane < lane_count; ++lane) {
            const std::size_t j = j_begin + lane;
            diag_buf[lane] =
                (j - 1u >= prev_begin && j - 1u <= prev_end) ? prev[j - 1u] : kInf;
            up_buf[lane] =
                (j >= prev_begin && j <= prev_end) ? prev[j] : kInf;
            ref_buf[lane] = ref[j - 1u];
        }
        for (std::size_t lane = lane_count; lane < kMaxVectorLanes; ++lane) {
            diag_buf[lane] = kInf;
            up_buf[lane] = kInf;
            ref_buf[lane] = 0;
        }

        const __mmask64 eq_mask = _mm512_cmpeq_epi8_mask(
            _mm512_set1_epi8(static_cast<char>(read[i - 1u])),
            _mm512_load_si512(reinterpret_cast<const __m512i*>(ref_buf.data())));
        const __m512i diag = _mm512_load_si512(reinterpret_cast<const __m512i*>(diag_buf.data()));
        const __m512i up = _mm512_load_si512(reinterpret_cast<const __m512i*>(up_buf.data()));
        const __m512i sub_cost = _mm512_mask_set1_epi32(_mm512_set1_epi32(1), static_cast<__mmask16>(eq_mask), 0);
        const __m512i best_diag = _mm512_add_epi32(diag, sub_cost);
        const __m512i best_up = _mm512_add_epi32(up, _mm512_set1_epi32(1));
        const __m512i base = _mm512_min_epi32(best_diag, best_up);
        _mm512_store_si512(reinterpret_cast<__m512i*>(base_buf.data()), base);

        int row_best = kInf;
        int left = row_has_zero ? static_cast<int>(i) : kInf;
        for (std::size_t lane = 0; lane < lane_count; ++lane) {
            int best = std::min(base_buf[lane], left + 1);
            curr[j_begin + lane] = best;
            left = best;
            row_best = std::min(row_best, best);
        }
        if (row_best > max_errors) {
            return row_best;
        }
        prev.swap(curr);
        prev_begin = curr_begin;
        prev_end = j_end;
    }

    return prev[m];
}

__attribute__((target("avx512f,avx512bw,avx512vl")))
int myers_avx512(const MyersQuery& query, std::string_view ref, int max_errors) {
    if (query.read_length == 0u) {
        return static_cast<int>(ref.size());
    }
    if (query.read_length > kMaxReadForMyers || max_errors > 7 || max_errors < 6) {
        return myers_popcnt(query, ref, max_errors);
    }
    std::vector<int> prev;
    std::vector<int> curr;
    return banded_score_only_avx512_impl(std::string_view(query.read_bases, query.read_length),
                                         ref,
                                         max_errors,
                                         prev,
                                         curr);
}

#endif

}  // namespace

MyersDispatch resolve_myers_dispatch() {
    MyersDispatch dispatch{};
    dispatch.selected = &myers_generic;
    dispatch.name = "generic";
#if defined(__x86_64__) || defined(__i386__)
    const SimdFeatures& features = detect_simd_features();
    // For mapper-speed's current public CLI (k<=3), the AVX512 banded verifier
    // is slower than the popcnt/bmi2 path because the DP band is too narrow to
    // amortize AVX512 setup and frequency costs. Keep the AVX512 kernel
    // available for future wider-band use, but prefer popcnt+bmi2 here.
    if (features.level == SimdLevel::kAvx2 && features.bmi2) {
        dispatch.selected = &myers_popcnt;
        dispatch.name = "popcnt+bmi2";
    } else if (features.level == SimdLevel::kAvx512 && features.bmi2) {
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
        query.read_bases[i] = read[i];
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
#if defined(__x86_64__) || defined(__i386__)
    const SimdFeatures& features = detect_simd_features();
    if (features.level == SimdLevel::kAvx512 &&
        features.avx512f && features.avx512bw && features.avx512vl &&
        static_cast<std::size_t>(max_errors * 2 + 1) >= 12u &&
        read.size() >= 96u && ref.size() >= 96u) {
        return banded_score_only_avx512_impl(read, ref, max_errors, prev, curr);
    }
#endif
    return banded_score_only_scalar_impl(read, ref, max_errors, prev, curr);
}

int banded_score_only(std::string_view read,
                      std::string_view ref,
                      int max_errors) {
    std::vector<int> prev;
    std::vector<int> curr;
    return banded_score_only(read, ref, max_errors, prev, curr);
}

}  // namespace mapper_speed

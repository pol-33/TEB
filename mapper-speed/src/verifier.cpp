#include "verifier.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>

namespace mapper_speed {

namespace {

struct PatternMasks128 {
    uint64_t lo[4]{};
    uint64_t hi[4]{};
};

PatternMasks128 build_pattern_masks(std::string_view read) {
    PatternMasks128 masks{};
    for (std::size_t i = 0; i < read.size(); ++i) {
        uint64_t bit = uint64_t{1} << (i & 63u);
        const std::size_t block = i >> 6u;
        uint64_t* dst = (block == 0u) ? masks.lo : masks.hi;
        switch (read[i]) {
            case 'A': dst[0] |= bit; break;
            case 'C': dst[1] |= bit; break;
            case 'G': dst[2] |= bit; break;
            case 'T': dst[3] |= bit; break;
            default: break;
        }
    }
    return masks;
}

template <bool UseBuiltinCpuHints>
int myers128_impl(std::string_view read, std::string_view ref, int max_errors) {
    if (read.empty()) {
        return static_cast<int>(ref.size());
    }
    if (read.size() > 128u) {
        return max_errors + 1;
    }

    const PatternMasks128 masks = build_pattern_masks(read);
    const bool has_hi = read.size() > 64u;

    uint64_t pv0 = ~uint64_t{0};
    uint64_t mv0 = 0;
    uint64_t pv1 = has_hi ? ~uint64_t{0} : 0;
    uint64_t mv1 = 0;
    const uint64_t last_bit0 = uint64_t{1} << ((std::min<std::size_t>(read.size(), 64u) - 1u) & 63u);
    const uint64_t last_bit1 = has_hi ? (uint64_t{1} << ((read.size() - 65u) & 63u)) : 0;
    int score = static_cast<int>(read.size());

    for (char base : ref) {
        const uint64_t* block0 = masks.lo;
        const uint64_t* block1 = masks.hi;
        std::size_t idx = 0;
        switch (base) {
            case 'A': idx = 0; break;
            case 'C': idx = 1; break;
            case 'G': idx = 2; break;
            case 'T': idx = 3; break;
            default: idx = 4; break;
        }

        const uint64_t eq0 = (idx < 4u) ? block0[idx] : 0;
        const uint64_t xv0 = eq0 | mv0;
        const uint64_t xh0 = (((eq0 & pv0) + pv0) ^ pv0) | eq0;
        uint64_t ph0 = mv0 | ~(xh0 | pv0);
        uint64_t mh0 = pv0 & xh0;

        uint64_t carry_ph = (ph0 >> 63u) & 1u;
        uint64_t carry_mh = (mh0 >> 63u) & 1u;
        ph0 = (ph0 << 1u) | 1u;
        mh0 <<= 1u;
        pv0 = mh0 | ~(xv0 | ph0);
        mv0 = ph0 & xv0;

        if (has_hi) {
            const uint64_t eq1 = (idx < 4u) ? block1[idx] : 0;
            const uint64_t xv1 = eq1 | mv1;
            const uint64_t xh1 = (((eq1 & pv1) + pv1 + carry_ph) ^ pv1) | eq1;
            uint64_t ph1 = mv1 | ~(xh1 | pv1);
            uint64_t mh1 = pv1 & xh1;
            ph1 = (ph1 << 1u) | carry_ph;
            mh1 = (mh1 << 1u) | carry_mh;
            pv1 = mh1 | ~(xv1 | ph1);
            mv1 = ph1 & xv1;
            if ((ph1 & last_bit1) != 0u) {
                ++score;
            } else if ((mh1 & last_bit1) != 0u) {
                --score;
            }
        } else {
            if ((ph0 & last_bit0) != 0u) {
                ++score;
            } else if ((mh0 & last_bit0) != 0u) {
                --score;
            }
        }

        if (UseBuiltinCpuHints && score > max_errors + static_cast<int>(ref.size())) {
            return score;
        }
    }

    return score;
}

int myers_generic(std::string_view read, std::string_view ref, int max_errors) {
    return myers128_impl<false>(read, ref, max_errors);
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("popcnt,bmi2")))
int myers_popcnt(std::string_view read, std::string_view ref, int max_errors) {
    return myers128_impl<true>(read, ref, max_errors);
}

__attribute__((target("avx512bw,avx512vl,avx512vbmi2,popcnt,bmi2")))
int myers_avx512(std::string_view read, std::string_view ref, int max_errors) {
    return myers128_impl<true>(read, ref, max_errors);
}
#endif

}  // namespace

MyersDispatch resolve_myers_dispatch() {
    MyersDispatch dispatch{};
    dispatch.selected = &myers_generic;
    dispatch.name = "generic";
#if defined(__x86_64__) || defined(__i386__)
    if (__builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512vl") &&
        __builtin_cpu_supports("avx512vbmi2") && __builtin_cpu_supports("bmi2")) {
        dispatch.selected = &myers_avx512;
        dispatch.name = "avx512";
    } else if (__builtin_cpu_supports("bmi2") && __builtin_cpu_supports("popcnt")) {
        dispatch.selected = &myers_popcnt;
        dispatch.name = "popcnt+bmi2";
    }
#endif
    return dispatch;
}

int bounded_edit_distance(const MyersDispatch& dispatch,
                          std::string_view read,
                          std::string_view ref,
                          int max_errors) {
    return dispatch.selected(read, ref, max_errors);
}

}  // namespace mapper_speed

#include "simd_dispatch.hpp"

#if !MAPPER_MEMORY_HAVE_AVX512_IMPL
#error "simd_avx512.cpp requires MAPPER_MEMORY_HAVE_AVX512_IMPL=1"
#endif

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>

#include <immintrin.h>

namespace mapper_memory::simd {

namespace {

using NibbleLutStorage = std::array<std::array<uint8_t, 64>, 4>;

const NibbleLutStorage& nibble_luts() {
    static const NibbleLutStorage values = [] {
        NibbleLutStorage tables{};
        for (std::size_t rank = 0; rank < tables.size(); ++rank) {
            std::array<uint8_t, 16> nibble_counts{};
            for (std::size_t nibble = 0; nibble < nibble_counts.size(); ++nibble) {
                uint8_t count = 0;
                count += static_cast<uint8_t>((nibble & 0x3U) == rank);
                count += static_cast<uint8_t>(((nibble >> 2U) & 0x3U) == rank);
                nibble_counts[nibble] = count;
            }
            for (std::size_t lane = 0; lane < 4; ++lane) {
                for (std::size_t nibble = 0; nibble < nibble_counts.size(); ++nibble) {
                    tables[rank][lane * 16U + nibble] = nibble_counts[nibble];
                }
            }
        }
        return tables;
    }();
    return values;
}

uint32_t count_partial_byte(uint8_t packed_rank,
                            uint8_t byte,
                            uint8_t first_symbol,
                            uint8_t end_symbol) {
    uint32_t count = 0;
    for (uint8_t symbol = first_symbol; symbol < end_symbol; ++symbol) {
        count += static_cast<uint32_t>(((byte >> (symbol * 2U)) & 0x3U) == packed_rank);
    }
    return count;
}

uint32_t reduce_counts(__m512i values) {
    alignas(64) uint64_t lanes[8];
    const __m512i sums = _mm512_sad_epu8(values, _mm512_setzero_si512());
    _mm512_store_si512(lanes, sums);

    uint32_t total = 0U;
    for (uint64_t lane : lanes) {
        total += static_cast<uint32_t>(lane);
    }
    return total;
}

bool range_contains(uint64_t start, uint64_t pos, uint64_t value) {
    return start <= value && value < pos;
}

}  // namespace

uint32_t count_packed_range_avx512(uint8_t packed_rank,
                                   const uint8_t* packed,
                                   uint64_t start,
                                   uint64_t pos,
                                   uint64_t primary_index) {
    if (packed == nullptr || start >= pos) {
        return 0U;
    }
    if (packed_rank > 3U) {
        return 0U;
    }

    uint64_t current = start;
    uint32_t count = 0U;

    if ((current & 3ULL) != 0ULL) {
        const uint64_t byte_index = current / 4ULL;
        const uint8_t first_symbol = static_cast<uint8_t>(current & 3ULL);
        const uint8_t end_symbol = static_cast<uint8_t>(std::min<uint64_t>(4ULL, pos - byte_index * 4ULL));
        count += count_partial_byte(packed_rank, packed[byte_index], first_symbol, end_symbol);
        current = std::min(pos, (byte_index + 1ULL) * 4ULL);
    }

    const __m512i low_mask = _mm512_set1_epi8(0x0F);
    const __m512i lut = _mm512_load_si512(reinterpret_cast<const __m512i*>(nibble_luts()[packed_rank].data()));

    uint64_t byte_index = current / 4ULL;
    const uint64_t full_byte_end = pos / 4ULL;
    for (; byte_index + 64ULL <= full_byte_end; byte_index += 64ULL) {
        const __m512i bytes = _mm512_loadu_si512(reinterpret_cast<const void*>(packed + byte_index));
        const __m512i low = _mm512_and_si512(bytes, low_mask);
        const __m512i high = _mm512_and_si512(_mm512_srli_epi16(bytes, 4), low_mask);
        const __m512i low_counts = _mm512_shuffle_epi8(lut, low);
        const __m512i high_counts = _mm512_shuffle_epi8(lut, high);
        count += reduce_counts(_mm512_add_epi8(low_counts, high_counts));
    }

    for (; byte_index < full_byte_end; ++byte_index) {
        const uint8_t byte = packed[byte_index];
        count += count_partial_byte(packed_rank, byte, 0U, 4U);
    }

    current = std::max(current, full_byte_end * 4ULL);
    if (current < pos) {
        const uint64_t tail_byte_index = current / 4ULL;
        const uint8_t end_symbol = static_cast<uint8_t>(pos - current);
        count += count_partial_byte(packed_rank, packed[tail_byte_index], 0U, end_symbol);
    }

    if (packed_rank == 0U && range_contains(start, pos, primary_index)) {
        --count;
    }
    return count;
}

}  // namespace mapper_memory::simd

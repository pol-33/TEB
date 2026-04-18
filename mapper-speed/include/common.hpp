#ifndef MAPPER_SPEED_COMMON_HPP
#define MAPPER_SPEED_COMMON_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <unistd.h>

#if defined(__x86_64__) || defined(__i386__)
#include <immintrin.h>
#endif

namespace mapper_speed {

constexpr uint32_t kSeedLength = 16;
constexpr uint32_t kOffsetPageShift = 12;
constexpr uint32_t kOffsetPageSize = 1u << kOffsetPageShift;
constexpr uint64_t kOffsetEntryCount = (uint64_t{1} << (2 * kSeedLength)) + 1u;
constexpr uint32_t kHighFreqSkip = 64;
constexpr uint32_t kHighFreqAllowFallback = 512;
constexpr uint32_t kMaxVerifyPerOrientation = 64;
constexpr uint32_t kCandidateTableSize = 8192;
constexpr uint32_t kMaxReadForMyers = 128;

inline uint64_t align_up(uint64_t value, uint64_t alignment) {
    return (value + alignment - 1u) & ~(alignment - 1u);
}

inline uint64_t fnv1a_extend(uint64_t hash, const void* data, std::size_t len) {
    const auto* bytes = static_cast<const uint8_t*>(data);
    for (std::size_t i = 0; i < len; ++i) {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
    return hash;
}

inline uint64_t fnv1a_init() {
    return 1469598103934665603ULL;
}

// Bulk byte mismatch counter using 64-bit SWAR.
// For each 8-byte chunk, XOR reveals differing bytes; zero-byte detection then
// counts equal bytes without branching, and popcount turns that into totals.
inline uint32_t count_byte_mismatches_swar(const char* lhs, const char* rhs, std::size_t len) {
    constexpr uint64_t kOnes = UINT64_C(0x0101010101010101);
    constexpr uint64_t kHighs = UINT64_C(0x8080808080808080);

    uint32_t mismatches = 0;
    std::size_t i = 0;
    for (; i + 8u <= len; i += 8u) {
        uint64_t a = 0;
        uint64_t b = 0;
        std::memcpy(&a, lhs + i, sizeof(uint64_t));
        std::memcpy(&b, rhs + i, sizeof(uint64_t));
        const uint64_t x = a ^ b;
        const uint64_t zero_high_bits = (x - kOnes) & ~x & kHighs;
        const uint32_t equal_bytes = static_cast<uint32_t>(__builtin_popcountll(zero_high_bits));
        mismatches += 8u - equal_bytes;
    }
    for (; i < len; ++i) {
        mismatches += static_cast<uint32_t>(lhs[i] != rhs[i]);
    }
    return mismatches;
}

#if defined(__x86_64__) || defined(__i386__)
__attribute__((target("avx2,popcnt")))
static inline uint32_t count_byte_mismatches_avx2(const char* lhs, const char* rhs, std::size_t len) {
    uint32_t mismatches = 0;
    std::size_t i = 0;
    for (; i + 32u <= len; i += 32u) {
        const __m256i a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(lhs + i));
        const __m256i b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(rhs + i));
        const __m256i eq = _mm256_cmpeq_epi8(a, b);
        const uint32_t mask = static_cast<uint32_t>(_mm256_movemask_epi8(eq));
        mismatches += 32u - static_cast<uint32_t>(__builtin_popcount(mask));
    }
    return mismatches + count_byte_mismatches_swar(lhs + i, rhs + i, len - i);
}

__attribute__((target("avx512bw,avx512vl,popcnt")))
static inline uint32_t count_byte_mismatches_avx512(const char* lhs, const char* rhs, std::size_t len) {
    uint32_t mismatches = 0;
    std::size_t i = 0;
    for (; i + 64u <= len; i += 64u) {
        const __m512i a = _mm512_loadu_si512(reinterpret_cast<const void*>(lhs + i));
        const __m512i b = _mm512_loadu_si512(reinterpret_cast<const void*>(rhs + i));
        const uint64_t mask = _mm512_cmpeq_epi8_mask(a, b);
        mismatches += 64u - static_cast<uint32_t>(__builtin_popcountll(mask));
    }
    return mismatches + count_byte_mismatches_swar(lhs + i, rhs + i, len - i);
}
#endif

inline uint32_t count_byte_mismatches_bulk(const char* lhs, const char* rhs, std::size_t len) {
#if defined(__x86_64__) || defined(__i386__)
    static const int level = []() {
        if (__builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512vl")) {
            return 2;
        }
        if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("popcnt")) {
            return 1;
        }
        return 0;
    }();
    if (level == 2) {
        return count_byte_mismatches_avx512(lhs, rhs, len);
    }
    if (level == 1) {
        return count_byte_mismatches_avx2(lhs, rhs, len);
    }
#endif
    return count_byte_mismatches_swar(lhs, rhs, len);
}

// Bulk GC counter: processes 8 bases per CPU cycle using 64-bit word tricks.
// For each byte, G(0x47) and C(0x43) have bits[1:0]==0b11; A, T, N do not.
// Strategy: extract bit0 and bit1 of every byte independently, AND them,
// then popcount — each byte with both bits set contributes exactly 1.
inline uint32_t count_gc_bases_bulk(const char* data, std::size_t len) {
    constexpr uint64_t kBit0 = UINT64_C(0x0101010101010101);
    uint32_t gc = 0;
    std::size_t i = 0;
    for (; i + 8u <= len; i += 8u) {
        uint64_t chunk = 0;
        std::memcpy(&chunk, data + i, sizeof(uint64_t));
        const uint64_t lo = chunk & kBit0;
        const uint64_t hi = (chunk >> 1u) & kBit0;
        gc += static_cast<uint32_t>(__builtin_popcountll(lo & hi));
    }
    for (; i < len; ++i) {
        gc += static_cast<uint32_t>((static_cast<uint8_t>(data[i]) & 0x3u) == 0x3u);
    }
    return gc;
}

inline double bytes_to_mebibytes(uint64_t bytes) {
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

inline uint64_t peak_rss_bytes() {
    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }
#if defined(__APPLE__)
    return static_cast<uint64_t>(usage.ru_maxrss);
#else
    return static_cast<uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

inline uint64_t physical_memory_bytes() {
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
    const long pages = sysconf(_SC_PHYS_PAGES);
    const long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0) {
        return static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_size);
    }
#endif
    return 0;
}

inline void throw_if(bool condition, const std::string& message) {
    if (condition) {
        throw std::runtime_error(message);
    }
}

}  // namespace mapper_speed

#endif

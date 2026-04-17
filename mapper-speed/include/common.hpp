#ifndef MAPPER_SPEED_COMMON_HPP
#define MAPPER_SPEED_COMMON_HPP

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sys/resource.h>

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

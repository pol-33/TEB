#ifndef MAPPER_MEMORY_MEMORY_STATS_HPP
#define MAPPER_MEMORY_MEMORY_STATS_HPP

#include <cstdint>
#include <string>

#if defined(__APPLE__)
#include <mach/mach.h>
#include <sys/resource.h>
#elif defined(__linux__)
#include <sys/resource.h>

#include <fstream>
#include <sstream>
#endif

namespace mapper_memory {

struct MemoryStats {
    uint64_t current_rss_bytes = 0;
    uint64_t peak_rss_bytes = 0;
};

inline MemoryStats read_memory_stats() {
    MemoryStats stats;

#if defined(__APPLE__)
    mach_task_basic_info_data_t info{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS) {
        stats.current_rss_bytes = static_cast<uint64_t>(info.resident_size);
    }

    struct rusage usage {};
    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        stats.peak_rss_bytes = static_cast<uint64_t>(usage.ru_maxrss);
    }
#elif defined(__linux__)
    std::ifstream status("/proc/self/status");
    std::string line;
    while (std::getline(status, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            std::istringstream iss(line.substr(6));
            uint64_t kb = 0;
            iss >> kb;
            stats.current_rss_bytes = kb * 1024ULL;
        } else if (line.rfind("VmPeak:", 0) == 0) {
            std::istringstream iss(line.substr(7));
            uint64_t kb = 0;
            iss >> kb;
            stats.peak_rss_bytes = kb * 1024ULL;
        }
    }

    if (stats.peak_rss_bytes == 0) {
        struct rusage usage {};
        if (getrusage(RUSAGE_SELF, &usage) == 0) {
            stats.peak_rss_bytes = static_cast<uint64_t>(usage.ru_maxrss) * 1024ULL;
        }
    }
#endif

    return stats;
}

inline double bytes_to_mebibytes(uint64_t bytes) {
    return static_cast<double>(bytes) / (1024.0 * 1024.0);
}

}  // namespace mapper_memory

#endif

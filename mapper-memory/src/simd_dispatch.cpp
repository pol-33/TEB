#include "simd_dispatch.hpp"

#include <array>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <stdexcept>
#include <string>

#if defined(__x86_64__) || defined(__i386__)
#include <cpuid.h>
#endif

namespace mapper_memory::simd {

namespace {

using CountFn = uint32_t (*)(uint8_t, const uint8_t*, uint64_t, uint64_t, uint64_t);

struct DispatchState {
    Backend backend = Backend::Scalar;
    CountFn count_fn = &count_packed_range_scalar;
    bool avx512_compiled = MAPPER_MEMORY_HAVE_AVX512_IMPL != 0;
    bool avx512_supported = false;
};

using ByteCountTable = std::array<std::array<uint8_t, 256>, 4>;

const ByteCountTable& byte_count_table() {
    static const ByteCountTable table = [] {
        ByteCountTable values{};
        for (std::size_t rank = 0; rank < values.size(); ++rank) {
            for (std::size_t byte = 0; byte < values[rank].size(); ++byte) {
                uint8_t count = 0;
                for (uint8_t symbol = 0; symbol < 4; ++symbol) {
                    count += static_cast<uint8_t>(((byte >> (symbol * 2U)) & 0x3U) == rank);
                }
                values[rank][byte] = count;
            }
        }
        return values;
    }();
    return table;
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

bool range_contains(uint64_t start, uint64_t pos, uint64_t value) {
    return start <= value && value < pos;
}

enum class RequestedBackend {
    Auto,
    Off,
    Avx512,
};

RequestedBackend requested_backend_from_env() {
    const char* env = std::getenv("MAPPER_SIMD");
    if (env == nullptr || *env == '\0' || std::strcmp(env, "auto") == 0) {
        return RequestedBackend::Auto;
    }
    if (std::strcmp(env, "off") == 0) {
        return RequestedBackend::Off;
    }
    if (std::strcmp(env, "avx512") == 0) {
        return RequestedBackend::Avx512;
    }
    throw std::runtime_error("invalid MAPPER_SIMD value: expected auto, off, or avx512");
}

#if defined(__x86_64__) || defined(__i386__)
uint64_t read_xcr0() {
    uint32_t eax = 0;
    uint32_t edx = 0;
    __asm__ volatile(".byte 0x0f, 0x01, 0xd0" : "=a"(eax), "=d"(edx) : "c"(0));
    return (static_cast<uint64_t>(edx) << 32U) | eax;
}

bool detect_avx512_support_on_host() {
    unsigned int eax = 0;
    unsigned int ebx = 0;
    unsigned int ecx = 0;
    unsigned int edx = 0;

    if (__get_cpuid(1, &eax, &ebx, &ecx, &edx) == 0) {
        return false;
    }
    if ((ecx & bit_OSXSAVE) == 0U) {
        return false;
    }

    constexpr uint64_t kRequiredXcr0Mask =
        (1ULL << 1U) | (1ULL << 2U) | (1ULL << 5U) | (1ULL << 6U) | (1ULL << 7U);
    if ((read_xcr0() & kRequiredXcr0Mask) != kRequiredXcr0Mask) {
        return false;
    }

    if (__get_cpuid_count(7, 0, &eax, &ebx, &ecx, &edx) == 0) {
        return false;
    }
    return (ebx & bit_AVX512F) != 0U && (ebx & bit_AVX512BW) != 0U;
}
#else
bool detect_avx512_support_on_host() {
    return false;
}
#endif

DispatchState make_dispatch_state() {
    DispatchState state;
    state.avx512_supported = avx512_supported_on_host();

    switch (requested_backend_from_env()) {
        case RequestedBackend::Off:
            return state;
        case RequestedBackend::Auto:
#if MAPPER_MEMORY_HAVE_AVX512_IMPL
            if (state.avx512_supported) {
                state.backend = Backend::Avx512;
                state.count_fn = &count_packed_range_avx512;
            }
#endif
            return state;
        case RequestedBackend::Avx512:
#if MAPPER_MEMORY_HAVE_AVX512_IMPL
            if (!state.avx512_supported) {
                throw std::runtime_error("MAPPER_SIMD=avx512 requested but AVX-512 is unavailable on this host");
            }
            state.backend = Backend::Avx512;
            state.count_fn = &count_packed_range_avx512;
            return state;
#else
            throw std::runtime_error("MAPPER_SIMD=avx512 requested but this binary was built without AVX-512 support");
#endif
    }

    return state;
}

const DispatchState& dispatch_state() {
    static std::once_flag once;
    static DispatchState state;
    std::call_once(once, [] { state = make_dispatch_state(); });
    return state;
}

}  // namespace

const char* backend_name(Backend backend) {
    switch (backend) {
        case Backend::Scalar:
            return "scalar";
        case Backend::Avx512:
            return "avx512";
    }
    return "unknown";
}

DispatchInfo resolved_dispatch() {
    const DispatchState& state = dispatch_state();
    return DispatchInfo{state.backend, state.avx512_compiled, state.avx512_supported};
}

uint32_t count_packed_range(uint8_t packed_rank,
                            const uint8_t* packed,
                            uint64_t start,
                            uint64_t pos,
                            uint64_t primary_index) {
    return dispatch_state().count_fn(packed_rank, packed, start, pos, primary_index);
}

uint32_t count_packed_range_scalar(uint8_t packed_rank,
                                   const uint8_t* packed,
                                   uint64_t start,
                                   uint64_t pos,
                                   uint64_t primary_index) {
    if (packed == nullptr || start >= pos) {
        return 0U;
    }
    if (packed_rank > 3U) {
        throw std::runtime_error("invalid packed rank for count_packed_range_scalar");
    }

    const ByteCountTable& table = byte_count_table();

    uint64_t current = start;
    uint32_t count = 0U;

    if ((current & 3ULL) != 0ULL) {
        const uint64_t byte_index = current / 4ULL;
        const uint8_t first_symbol = static_cast<uint8_t>(current & 3ULL);
        const uint8_t end_symbol = static_cast<uint8_t>(std::min<uint64_t>(4ULL, pos - byte_index * 4ULL));
        count += count_partial_byte(packed_rank, packed[byte_index], first_symbol, end_symbol);
        current = std::min<uint64_t>(pos, (byte_index + 1ULL) * 4ULL);
    }

    const uint64_t full_byte_end = pos / 4ULL;
    for (uint64_t byte_index = current / 4ULL; byte_index < full_byte_end; ++byte_index) {
        count += table[packed_rank][packed[byte_index]];
    }
    current = std::max<uint64_t>(current, full_byte_end * 4ULL);

    if (current < pos) {
        const uint64_t byte_index = current / 4ULL;
        const uint8_t end_symbol = static_cast<uint8_t>(pos - current);
        count += count_partial_byte(packed_rank, packed[byte_index], 0U, end_symbol);
    }

    if (packed_rank == 0U && range_contains(start, pos, primary_index)) {
        --count;
    }
    return count;
}

bool avx512_impl_available() {
    return MAPPER_MEMORY_HAVE_AVX512_IMPL != 0;
}

bool avx512_supported_on_host() {
    static const bool supported = detect_avx512_support_on_host();
    return supported;
}

}  // namespace mapper_memory::simd

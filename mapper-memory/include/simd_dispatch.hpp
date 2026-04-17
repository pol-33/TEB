#ifndef MAPPER_MEMORY_SIMD_DISPATCH_HPP
#define MAPPER_MEMORY_SIMD_DISPATCH_HPP

#include <cstdint>

#ifndef MAPPER_MEMORY_HAVE_AVX512_IMPL
#define MAPPER_MEMORY_HAVE_AVX512_IMPL 0
#endif

namespace mapper_memory::simd {

enum class Backend {
    Scalar,
    Avx512,
};

struct DispatchInfo {
    Backend backend = Backend::Scalar;
    bool avx512_compiled = false;
    bool avx512_supported = false;
};

const char* backend_name(Backend backend);
DispatchInfo resolved_dispatch();

uint32_t count_packed_range(uint8_t packed_rank,
                            const uint8_t* packed,
                            uint64_t start,
                            uint64_t pos,
                            uint64_t primary_index);

uint32_t count_packed_range_scalar(uint8_t packed_rank,
                                   const uint8_t* packed,
                                   uint64_t start,
                                   uint64_t pos,
                                   uint64_t primary_index);

bool avx512_impl_available();
bool avx512_supported_on_host();

#if MAPPER_MEMORY_HAVE_AVX512_IMPL
uint32_t count_packed_range_avx512(uint8_t packed_rank,
                                   const uint8_t* packed,
                                   uint64_t start,
                                   uint64_t pos,
                                   uint64_t primary_index);
#endif

}  // namespace mapper_memory::simd

#endif

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "fm_index.hpp"
#include "nucleotide.hpp"
#include "simd_dispatch.hpp"

namespace {

using mapper_memory::FMIndexView;

struct SyntheticChromosome {
    std::vector<uint8_t> packed_bwt;
    std::vector<uint32_t> sampled_occ;
    uint64_t text_length = 0;
    uint64_t primary_index = 0;
    uint32_t occ_sample = 1;
};

uint8_t dna_rank_from_packed(uint8_t packed_rank) {
    return static_cast<uint8_t>(mapper_memory::kARank + packed_rank);
}

SyntheticChromosome make_synthetic_chromosome(const std::vector<uint8_t>& packed_symbols,
                                              uint64_t primary_index,
                                              uint32_t occ_sample) {
    if (packed_symbols.empty()) {
        throw std::runtime_error("synthetic chromosome must be non-empty");
    }
    SyntheticChromosome chrom;
    chrom.text_length = static_cast<uint64_t>(packed_symbols.size());
    chrom.primary_index = primary_index;
    chrom.occ_sample = occ_sample;
    chrom.packed_bwt.assign(mapper_memory::packed_byte_count(chrom.text_length), 0U);

    for (uint64_t i = 0; i < chrom.text_length; ++i) {
        const uint8_t stored_rank = (i == primary_index)
            ? mapper_memory::kARank
            : dna_rank_from_packed(packed_symbols[static_cast<std::size_t>(i)]);
        mapper_memory::set_packed_rank(chrom.packed_bwt, i, stored_rank);
    }

    const std::size_t blocks = static_cast<std::size_t>(chrom.text_length / occ_sample) + 1U;
    chrom.sampled_occ.assign(blocks * 4U, 0U);
    std::array<uint32_t, 4> running{};

    for (uint64_t row = 0; row < chrom.text_length; ++row) {
        if (row != primary_index) {
            ++running[packed_symbols[static_cast<std::size_t>(row)]];
        }
        if ((row + 1ULL) % occ_sample == 0ULL) {
            const std::size_t block = static_cast<std::size_t>((row + 1ULL) / occ_sample);
            for (std::size_t rank = 0; rank < running.size(); ++rank) {
                chrom.sampled_occ[block * 4U + rank] = running[rank];
            }
        }
    }

    return chrom;
}

FMIndexView::ChromosomeView make_view(const SyntheticChromosome& chrom) {
    FMIndexView::ChromosomeView view;
    view.name = "synthetic";
    view.length = chrom.text_length;
    view.text_length = chrom.text_length;
    view.primary_index = chrom.primary_index;
    view.occ_sample = chrom.occ_sample;
    view.sa_sample = 1;
    view.packed_bwt = chrom.packed_bwt.data();
    view.packed_bwt_bytes = chrom.packed_bwt.size();
    view.sampled_occ = chrom.sampled_occ.data();
    view.sampled_occ_entries = chrom.sampled_occ.size();
    return view;
}

uint32_t naive_occ(const std::vector<uint8_t>& packed_symbols,
                   uint64_t primary_index,
                   uint8_t rank,
                   uint64_t pos) {
    if (rank == mapper_memory::kSentinelRank) {
        return static_cast<uint32_t>(pos > primary_index ? 1U : 0U);
    }
    if (rank == mapper_memory::kSeparatorRank || rank > mapper_memory::kTRank) {
        return 0U;
    }
    uint32_t count = 0U;
    for (uint64_t i = 0; i < pos && i < packed_symbols.size(); ++i) {
        if (i == primary_index) {
            continue;
        }
        count += static_cast<uint32_t>(dna_rank_from_packed(packed_symbols[static_cast<std::size_t>(i)]) == rank);
    }
    return count;
}

void assert_occ_matches_naive(const std::vector<uint8_t>& packed_symbols,
                              uint64_t primary_index,
                              uint32_t occ_sample) {
    const SyntheticChromosome chrom = make_synthetic_chromosome(packed_symbols, primary_index, occ_sample);
    const FMIndexView::ChromosomeView view = make_view(chrom);

    for (uint8_t rank = mapper_memory::kSentinelRank; rank <= mapper_memory::kTRank; ++rank) {
        for (uint64_t pos = 0; pos <= chrom.text_length; ++pos) {
            const uint32_t expected = naive_occ(packed_symbols, primary_index, rank, pos);
            const uint32_t actual = view.occ(rank, pos);
            if (actual != expected) {
                throw std::runtime_error(
                    "occ mismatch: rank=" + std::to_string(rank) +
                    " pos=" + std::to_string(pos) +
                    " expected=" + std::to_string(expected) +
                    " actual=" + std::to_string(actual));
            }
        }
    }
}

#if MAPPER_MEMORY_HAVE_AVX512_IMPL
void assert_backends_match(const std::vector<uint8_t>& packed_symbols,
                           uint64_t primary_index) {
    const SyntheticChromosome chrom = make_synthetic_chromosome(packed_symbols, primary_index, 1U);
    const uint64_t n = chrom.text_length;

    for (uint8_t packed_rank = 0; packed_rank < 4; ++packed_rank) {
        for (uint64_t start = 0; start <= n; ++start) {
            for (uint64_t pos = start; pos <= n; ++pos) {
                const uint32_t scalar = mapper_memory::simd::count_packed_range_scalar(
                    packed_rank, chrom.packed_bwt.data(), start, pos, primary_index);
#if MAPPER_MEMORY_HAVE_AVX512_IMPL
                const uint32_t avx512 = mapper_memory::simd::count_packed_range_avx512(
                    packed_rank, chrom.packed_bwt.data(), start, pos, primary_index);
                if (scalar != avx512) {
                    throw std::runtime_error(
                        "backend mismatch: rank=" + std::to_string(packed_rank) +
                        " start=" + std::to_string(start) +
                        " pos=" + std::to_string(pos) +
                        " scalar=" + std::to_string(scalar) +
                        " avx512=" + std::to_string(avx512));
                }
#else
                (void)scalar;
#endif
            }
        }
    }
}
#endif

void run_occ_regression_suite() {
    const std::vector<uint32_t> occ_samples = {1U, 3U, 4U, 5U, 8U, 16U, 32U, 256U};
    const std::vector<std::vector<uint8_t>> fixed_sequences = {
        {0U},
        {0U, 1U, 2U, 3U},
        {3U, 3U, 2U, 1U, 0U, 0U, 1U},
        {0U, 1U, 2U, 3U, 0U, 1U, 2U, 3U, 0U},
        {3U, 2U, 1U, 0U, 3U, 2U, 1U, 0U, 3U, 2U, 1U},
    };

    for (const auto& sequence : fixed_sequences) {
        for (uint64_t primary = 0; primary < sequence.size(); ++primary) {
            for (uint32_t occ_sample : occ_samples) {
                assert_occ_matches_naive(sequence, primary, occ_sample);
            }
        }
    }

    std::mt19937 rng(123456789U);
    std::uniform_int_distribution<int> rank_dist(0, 3);
    std::uniform_int_distribution<int> length_dist(1, 129);

    for (int iter = 0; iter < 48; ++iter) {
        const int length = length_dist(rng);
        std::vector<uint8_t> sequence(static_cast<std::size_t>(length), 0U);
        for (uint8_t& symbol : sequence) {
            symbol = static_cast<uint8_t>(rank_dist(rng));
        }

        for (uint64_t primary = 0; primary < sequence.size(); primary += std::max<uint64_t>(1ULL, sequence.size() / 5ULL)) {
            for (uint32_t occ_sample : occ_samples) {
                assert_occ_matches_naive(sequence, primary, occ_sample);
            }
            if (primary + 1U == sequence.size()) {
                break;
            }
        }
    }

    std::cout << "[simd-selftest] occ regression suite passed\n";
}

void run_backend_comparison_suite() {
#if MAPPER_MEMORY_HAVE_AVX512_IMPL
    if (!mapper_memory::simd::avx512_supported_on_host()) {
        std::cout << "[simd-selftest] AVX-512 compiled but unavailable on this host; skipping backend comparison\n";
        return;
    }

    std::mt19937 rng(424242U);
    std::uniform_int_distribution<int> rank_dist(0, 3);

    for (std::size_t length : {1U, 2U, 3U, 4U, 5U, 7U, 15U, 31U, 63U, 64U, 65U, 127U}) {
        std::vector<uint8_t> sequence(length, 0U);
        for (uint8_t& symbol : sequence) {
            symbol = static_cast<uint8_t>(rank_dist(rng));
        }
        for (uint64_t primary = 0; primary < sequence.size(); primary += std::max<uint64_t>(1ULL, sequence.size() / 4ULL)) {
            assert_backends_match(sequence, primary);
            if (primary + 1U == sequence.size()) {
                break;
            }
        }
    }

    std::cout << "[simd-selftest] scalar vs AVX-512 backend comparison passed\n";
#else
    std::cout << "[simd-selftest] AVX-512 backend not built; backend comparison skipped\n";
#endif
}

}  // namespace

int main() {
    try {
        const mapper_memory::simd::DispatchInfo dispatch = mapper_memory::simd::resolved_dispatch();
        std::cout << "[simd-selftest] selected backend: "
                  << mapper_memory::simd::backend_name(dispatch.backend) << '\n';

        run_occ_regression_suite();
        run_backend_comparison_suite();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "simd_selftest error: " << e.what() << '\n';
        return 1;
    }
}

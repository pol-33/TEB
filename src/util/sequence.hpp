#pragma once

#include <cstdint>
#ifdef __AVX2__
#  include <immintrin.h>
#endif

// ---------------------------------------------------------------------------
// EXPAND_LUT[256]
//
// Maps one byte of n_mask (1 bit per base, LSB = base 0) to a uint16_t
// with 2 bits per base so the result aligns with the 2-bit packed encoding.
// Each input bit b is "doubled": bit b → bits [2b+1 : 2b] = 0b11 in output.
//
// Example:  0b00000001  →  0b0000000000000011   (base 0 is N)
//           0b00000010  →  0b0000000000001100   (base 1 is N)
//           0b10000001  →  0b1100000000000011   (bases 0 and 7 are N)
// ---------------------------------------------------------------------------
namespace detail_seq {

struct ExpandLutHolder {
    uint16_t v[256];
    constexpr ExpandLutHolder() noexcept : v() {
        for (int i = 0; i < 256; ++i) {
            uint16_t out = 0;
            for (int b = 0; b < 8; ++b)
                if (i & (1 << b))
                    out = static_cast<uint16_t>(out | (3u << (2 * b)));
            v[i] = out;
        }
    }
};

static constexpr ExpandLutHolder EXPAND_LUT_HOLDER{};

} // namespace detail_seq

// Public reference to the table.
static constexpr const uint16_t (&EXPAND_LUT)[256] = detail_seq::EXPAND_LUT_HOLDER.v;

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Extract 8 consecutive n_mask bits starting at absolute base position p.
// Correctly handles the case where the 8 bits cross a 64-bit word boundary
// (only possible when p % 64 > 56).
static inline uint8_t nmask_extract8(const uint64_t* nm, uint32_t p) noexcept {
    const uint32_t w   = p / 64;
    const uint32_t bit = p % 64;
    const uint64_t lo  = nm[w] >> bit;
    if (bit <= 56) return static_cast<uint8_t>(lo & 0xFFu);
    // bits_in_lo = 64 - bit < 8  →  borrow remaining bits from next word
    return static_cast<uint8_t>((lo | (nm[w + 1] << (64u - bit))) & 0xFFu);
}

// Extract a genome packed byte (4 consecutive bases) starting at absolute
// base position p.  Handles non-4-aligned p by merging two adjacent bytes.
static inline uint8_t genome_extract_byte(const uint8_t* pk, uint32_t p) noexcept {
    const uint32_t bi    = p / 4;
    const uint32_t shift = (p % 4) * 2;   // 0, 2, 4, or 6
    if (shift == 0) return pk[bi];
    return static_cast<uint8_t>((pk[bi] >> shift) | (pk[bi + 1] << (8u - shift)));
}

// ---------------------------------------------------------------------------
// count_mismatches
//
// Counts the number of base-level mismatches between
//   genome_packed[genome_pos .. genome_pos + read_len)
// and
//   read_packed[0 .. read_len)
//
// Positions where the genome has N (bit set in genome_n_mask) are treated as
// wildcards and never counted as mismatches.
//
// Encoding (shared by genome and read):
//   A = 0b00,  C = 0b01,  G = 0b10,  T = 0b11
//   4 bases per byte, LSB-first (base i is at bits [(i%4)*2 + 1 : (i%4)*2]
//   of byte i/4).
//   N is stored as A (0b00) in the packed array; the n_mask is the
//   authoritative record of N positions.
//
// AVX2 path  (#ifdef __AVX2__): 128 bases per iteration.
// Scalar fallback (always compiled):  4 bases per iteration.
// ---------------------------------------------------------------------------
inline int count_mismatches(
    const uint8_t*  genome_packed,
    const uint64_t* genome_n_mask,
    uint32_t        genome_pos,
    const uint8_t*  read_packed,
    uint32_t        read_len
) {
    int      mismatches = 0;
    uint32_t pos        = 0;   // current offset within the read (in bases)

#ifdef __AVX2__
    // -----------------------------------------------------------------------
    // AVX2 path – 128 bases (32 packed bytes) per iteration.
    //
    // For each 128-base window:
    //   1. Build 32-byte genome buffer (handles non-4-aligned genome_pos).
    //   2. Expand 16 n_mask bytes → 32 bytes using EXPAND_LUT.
    //   3. XOR genome and read packed bytes.
    //   4. Zero out N positions:  ~nmask AND xor.
    //   5. If non-zero, count non-zero 2-bit pairs = mismatches.
    // -----------------------------------------------------------------------
    while (pos + 128 <= read_len) {
        const uint32_t gs = genome_pos + pos;

        // --- genome window -----------------------------------------------
        alignas(32) uint8_t g_buf[32];
        for (int j = 0; j < 32; ++j)
            g_buf[j] = genome_extract_byte(genome_packed, gs + static_cast<uint32_t>(j) * 4);

        // --- expanded n_mask window (16 n_mask bytes → 32 packed-mask bytes)
        alignas(32) uint8_t nm_buf[32];
        for (int j = 0; j < 16; ++j) {
            const uint16_t e = EXPAND_LUT[nmask_extract8(genome_n_mask, gs + static_cast<uint32_t>(j) * 8)];
            nm_buf[2 * j]     = static_cast<uint8_t>(e & 0xFFu);
            nm_buf[2 * j + 1] = static_cast<uint8_t>(e >> 8);
        }

        __m256i gv  = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(g_buf));
        __m256i rv  = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(read_packed + pos / 4));
        __m256i xv  = _mm256_xor_si256(gv, rv);
        __m256i nmv = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(nm_buf));
        // ~nmv & xv  →  bits set only where genome != read AND genome is not N
        __m256i mv  = _mm256_andnot_si256(nmv, xv);

        if (!_mm256_testz_si256(mv, mv)) {
            // Count non-zero 2-bit pairs across all 256 bits.
            alignas(32) uint64_t lanes[4];
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(lanes), mv);
            for (int j = 0; j < 4; ++j) {
                uint64_t v = lanes[j];
                v = (v | (v >> 1)) & UINT64_C(0x5555555555555555);
                mismatches += __builtin_popcountll(v);
            }
        }

        pos += 128;
    }
#endif // __AVX2__

    // -----------------------------------------------------------------------
    // Scalar path – 4 bases (1 packed byte) per iteration.
    //
    //   XOR:      g_byte ^ r_byte   – bits differ where bases differ
    //   n_mask:   4-bit nibble from genome_n_mask
    //   EXPAND:   EXPAND_LUT[nibble] lower byte → 8-bit packed mask
    //   mask out: xor & ~nmask8     – zero out N positions
    //   count:    non-zero 2-bit pairs in the masked XOR byte
    // -----------------------------------------------------------------------
    while (pos + 4 <= read_len) {
        const uint32_t gs = genome_pos + pos;

        const uint8_t g_byte = genome_extract_byte(genome_packed, gs);
        const uint8_t r_byte = read_packed[pos / 4];
        const uint8_t xored  = g_byte ^ r_byte;

        // Extract 4-bit n_mask nibble for bases gs..gs+3.
        const uint32_t w   = gs / 64;
        const uint32_t b   = gs % 64;
        uint8_t nm4;
        if (b <= 60) {
            nm4 = static_cast<uint8_t>((genome_n_mask[w] >> b) & 0xFu);
        } else {
            // Cross-word boundary (b = 61, 62, or 63).
            nm4 = static_cast<uint8_t>(
                ((genome_n_mask[w] >> b) | (genome_n_mask[w + 1] << (64u - b))) & 0xFu);
        }

        // Lower byte of EXPAND_LUT[nm4] is the 8-bit packed mask for 4 bases.
        const uint8_t nm8   = static_cast<uint8_t>(EXPAND_LUT[nm4] & 0xFFu);
        const uint8_t masked = static_cast<uint8_t>(xored & ~nm8);

        // Fold each 2-bit pair to its lower bit, then popcount.
        const uint8_t v = static_cast<uint8_t>((masked | (masked >> 1u)) & 0x55u);
        mismatches += __builtin_popcount(v);

        pos += 4;
    }

    // -----------------------------------------------------------------------
    // Tail – up to 3 remaining bases, one at a time.
    // -----------------------------------------------------------------------
    while (pos < read_len) {
        const uint32_t gp = genome_pos + pos;
        const uint8_t  gb = (genome_packed[gp  / 4] >> ((gp  % 4) * 2)) & 3u;
        const uint8_t  rb = (read_packed[pos / 4]   >> ((pos % 4) * 2)) & 3u;
        const bool     gn = (genome_n_mask[gp / 64] >> (gp % 64)) & 1u;
        if (!gn && gb != rb) ++mismatches;
        ++pos;
    }

    return mismatches;
}

/*
 * ===========================================================================
 * SELF-TEST
 * ===========================================================================
 *
 * Build and run with:
 *   g++ -std=c++17 -O2 -mavx2 -DSEQUENCE_SELF_TEST \
 *       -I src/util src/util/sequence.hpp -o seq_test && ./seq_test
 *
 * (or without -mavx2 to exercise only the scalar path)
 * ===========================================================================
 *
 * #include "sequence.hpp"
 * #include <cassert>
 * #include <cstdio>
 * #include <cstring>
 *
 * // Pack a plain ASCII sequence into the 2-bit LSB-first format used by
 * // GenomeStorage.  N is encoded as A (00) and its position is recorded in
 * // n_mask.
 * static void pack(const char* seq, int len,
 *                  uint8_t* packed, uint64_t* n_mask) {
 *     const int n_words = (len + 63) / 64;
 *     memset(packed,  0, (len + 3) / 4);
 *     memset(n_mask,  0, n_words * 8);
 *     for (int i = 0; i < len; ++i) {
 *         char c = seq[i];
 *         uint8_t code = 0;
 *         if      (c=='C'||c=='c') code = 1;
 *         else if (c=='G'||c=='g') code = 2;
 *         else if (c=='T'||c=='t') code = 3;
 *         else if (c=='N'||c=='n') { n_mask[i/64] |= (uint64_t)1 << (i%64); }
 *         packed[i/4] |= (uint8_t)(code << ((i%4)*2));
 *     }
 * }
 *
 * // -----------------------------------------------------------------------
 * // Test 1 – short (8 bases), exercises scalar path only.
 * //
 * //   Genome (genome_pos=0): A C G T N A C G
 * //                           0 1 2 3 4 5 6 7
 * //   Read:                   T C G T T A C G
 * //
 * //   pos 0: A(00) vs T(11) → mismatch
 * //   pos 4: N vs T         → genome N → ignored
 * //   All others:           → match
 * //   Expected: 1 mismatch
 * //
 * // Without masking pos 4 would yield 2 mismatches (false positive at pos 4
 * // because N is stored as A=00 and read has T=11, XOR=11 ≠ 0).
 * // -----------------------------------------------------------------------
 * static void test_short() {
 *     const char genome_seq[] = "ACGTNACG";   // 8 bases, N at pos 4
 *     const char read_seq[]   = "TCGTTACG";   // mismatch at pos 0; T at pos 4
 *     const int  LEN = 8;
 *
 *     uint8_t  gp[2]  = {};
 *     uint64_t nm[1]  = {};
 *     uint8_t  rp[2]  = {};
 *     uint64_t dummy  = 0;  // read has no Ns
 *
 *     pack(genome_seq, LEN, gp, nm);
 *     pack(read_seq,   LEN, rp, &dummy);
 *
 *     int mm = count_mismatches(gp, nm, 0, rp, LEN);
 *     assert(mm == 1 && "test_short: expected 1 mismatch");
 *     printf("test_short PASS  (mm=%d, expected 1)\n", mm);
 * }
 *
 * // -----------------------------------------------------------------------
 * // Test 2 – genome_pos != 0 (non-4-aligned), scalar path.
 * //
 * //   Genome buffer (pos 0..11): A A A A A A C G T N A C
 * //   genome_pos = 6, so we compare genome[6..11] = C G T N A C
 * //   Read (6 bases):                                 C G T T A C
 * //
 * //   pos 0 (genome 6): C vs C → match
 * //   pos 1 (genome 7): G vs G → match
 * //   pos 2 (genome 8): T vs T → match
 * //   pos 3 (genome 9): N vs T → ignored (genome N)
 * //   pos 4 (genome10): A vs A → match
 * //   pos 5 (genome11): C vs C → match
 * //   Expected: 0 mismatches  (N at genome pos 9 is correctly masked)
 * // -----------------------------------------------------------------------
 * static void test_nonzero_genome_pos() {
 *     const char genome_seq[] = "AAAAAACGTNAC";  // 12 bases, N at pos 9
 *     const char read_seq[]   = "CGTTAC";         // 6 bases
 *     const int  GLEN = 12, RLEN = 6;
 *
 *     uint8_t  gp[3]  = {};
 *     uint64_t nm[1]  = {};
 *     uint8_t  rp[2]  = {};
 *     uint64_t dummy  = 0;
 *
 *     pack(genome_seq, GLEN, gp, nm);
 *     pack(read_seq,   RLEN, rp, &dummy);
 *
 *     int mm = count_mismatches(gp, nm, 6, rp, RLEN);
 *     assert(mm == 0 && "test_nonzero_genome_pos: expected 0 mismatches");
 *     printf("test_nonzero_genome_pos PASS  (mm=%d, expected 0)\n", mm);
 * }
 *
 * // -----------------------------------------------------------------------
 * // Test 3 – 256 bases, exercises AVX2 path (first 128-base chunk) and
 * //          the scalar tail for the remaining 128 bases.
 * //
 * //   Genome (pos 0..255): "ACGT" repeated 64 times, with N at positions
 * //                         60 and 190.
 * //   Read:                 identical copy except mismatches at positions
 * //                         5 (C→G), 130 (G→A), and 200 (T→A).
 * //                         Position 60 in read is T (would mismatch N, but
 * //                         should be ignored).
 * //
 * //   Mismatch at pos  5:  genome C(01) vs read G(10) → differ (+1)
 * //   Mismatch at pos 60:  genome N vs read T         → ignored (genome N)
 * //   Mismatch at pos130:  genome G(10) vs read A(00) → differ (+1)
 * //   Mismatch at pos190:  genome N vs read A         → ignored (genome N)
 *   Mismatch at pos200:  genome A(00) vs read T(11) → differ (+1)
 * //   Expected: 3 mismatches
 * //
 * //   The scalar-only count (computed by temporarily #undefing AVX2, or by
 * //   using a reference loop) must equal the AVX2+scalar count.
 * // -----------------------------------------------------------------------
 * static void test_long_avx2_vs_scalar() {
 *     const int LEN = 256;
 *     char genome_seq[LEN+1], read_seq[LEN+1];
 *     const char* pat = "ACGT";
 *     for (int i = 0; i < LEN; ++i) genome_seq[i] = pat[i % 4];
 *     memcpy(read_seq, genome_seq, LEN);
 *     genome_seq[LEN] = read_seq[LEN] = '\0';
 *
 *     // Introduce Ns in genome (stored as A=00 in packed, bit set in n_mask)
 *     genome_seq[60]  = 'N';
 *     genome_seq[190] = 'N';
 *
 *     // Introduce mismatches in read
 *     read_seq[5]   = 'G';   // C->G  (5%4=1 -> C, changed to G)
 *     read_seq[60]  = 'T';   // N position in genome -> ignored
 *     read_seq[130] = 'A';   // G->A  (130%4=2 -> G, changed to A)
 *     read_seq[190] = 'A';   // N position in genome -> ignored
 *     read_seq[200] = 'T';   // A->T  (200%4=0 -> A, changed to T)
 *
 *     uint8_t  gp[(LEN+3)/4];
 *     uint64_t nm[(LEN+63)/64 + 1];  // +1 for cross-word safety
 *     uint8_t  rp[(LEN+3)/4];
 *     uint64_t dummy[(LEN+63)/64 + 1];
 *     memset(gp, 0, sizeof(gp));
 *     memset(nm, 0, sizeof(nm));
 *     memset(rp, 0, sizeof(rp));
 *     memset(dummy, 0, sizeof(dummy));
 *
 *     pack(genome_seq, LEN, gp, nm);
 *     pack(read_seq,   LEN, rp, dummy);
 *
 *     int mm = count_mismatches(gp, nm, 0, rp, LEN);
 *     assert(mm == 3 && "test_long_avx2_vs_scalar: expected 3 mismatches");
 *     printf("test_long_avx2_vs_scalar PASS  (mm=%d, expected 3)\n", mm);
 *
 *     // Also verify scalar-only by computing with the reference loop.
 *     int scalar_mm = 0;
 *     for (int i = 0; i < LEN; ++i) {
 *         bool gn = (nm[i/64] >> (i%64)) & 1u;
 *         if (gn) continue;
 *         uint8_t gb = (gp[i/4] >> ((i%4)*2)) & 3u;
 *         uint8_t rb = (rp[i/4] >> ((i%4)*2)) & 3u;
 *         if (gb != rb) ++scalar_mm;
 *     }
 *     assert(scalar_mm == 3 && "scalar reference must also give 3");
 *     assert(mm == scalar_mm && "AVX2+scalar path must match scalar reference");
 *     printf("test_long_avx2_vs_scalar – AVX2 path matches scalar reference PASS\n");
 * }
 *
 * int main() {
 *     test_short();
 *     test_nonzero_genome_pos();
 *     test_long_avx2_vs_scalar();
 *     printf("All tests PASSED\n");
 *     return 0;
 * }
 *
 * ===========================================================================
 * END SELF-TEST
 * ===========================================================================
 */

// test_kmer_roundtrip.cpp
// Standalone roundtrip test: build → save → load(SPEED) → load(MEMORY) → query
//
// Compile (from project root):
//   g++ -std=c++17 -O3 -march=native \
//       -I src/index -I src/genome -I src/io -I src/util -I vendor \
//       src/index/kmer_index.cpp src/index/index.cpp \
//       src/genome/genome.cpp src/io/fasta.cpp \
//       test_kmer_roundtrip.cpp -o /tmp/test_kmer_roundtrip

#include "src/index/kmer_index.hpp"

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Write a minimal FASTA file with a single sequence.
static void write_fasta(const std::string& path, const std::string& seq) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot create " + path);
    f << ">test_seq\n" << seq << "\n";
}

// Pack `read` (ASCII, uppercase ACGT) into 2-bit format, LSB-first.
// Returns vector of bytes; caller must ensure read length is correct.
static uint8_t base_enc(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0;
    }
}

static std::vector<uint8_t> pack_read(const std::string& read) {
    const size_t bytes = (read.size() + 3) / 4;
    std::vector<uint8_t> out(bytes, 0);
    for (size_t i = 0; i < read.size(); ++i)
        out[i / 4] |= static_cast<uint8_t>(base_enc(read[i]) << ((i % 4) * 2));
    return out;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    constexpr size_t GENOME_LEN = 10000;
    constexpr size_t K          = 14;

    // -----------------------------------------------------------------------
    // Build a synthetic genome of length 10000.
    //
    // Layout:
    //   positions [   0 ..  13]: MARKER k-mer "ACGTACGTACGTAC"
    //   positions [1000 .. 1013]: MARKER k-mer repeated
    //   all other positions: 'G'
    //
    // The MARKER never appears in a pure-G background, so query returns
    // exactly {0, 1000}.
    // -----------------------------------------------------------------------
    const std::string MARKER = "ACGTACGTACGTAC";  // len == K == 14
    assert(MARKER.size() == K);

    std::string genome(GENOME_LEN, 'G');
    genome.replace(0,    K, MARKER);
    genome.replace(1000, K, MARKER);

    const std::string fasta_path = "/tmp/teb_test_genome.fna";
    const std::string idx_path   = "/tmp/teb_test.idx";
    write_fasta(fasta_path, genome);

    // -----------------------------------------------------------------------
    // build()
    // -----------------------------------------------------------------------
    std::cout << "=== build() ===\n";
    KmerIndex idx(K);
    idx.build(fasta_path);

    {
        const auto positions = idx.query(MARKER);
        std::vector<uint32_t> expected = {0, 1000};
        std::sort(const_cast<std::vector<uint32_t>&>(positions).begin(),
                  const_cast<std::vector<uint32_t>&>(positions).end());
        assert(positions == expected && "build: query mismatch");
        std::cout << "  query post-build: OK (positions: ";
        for (auto p : positions) std::cout << p << " ";
        std::cout << ")\n";
    }

    // -----------------------------------------------------------------------
    // save()
    // -----------------------------------------------------------------------
    std::cout << "=== save() ===\n";
    idx.save(idx_path);
    std::cout << "  saved to " << idx_path << "\n";

    // -----------------------------------------------------------------------
    // load(SPEED) + query + verify_match
    // -----------------------------------------------------------------------
    std::cout << "=== load(SPEED) ===\n";
    {
        KmerIndex idx2(K);
        idx2.load(idx_path, LoadMode::SPEED);

        auto positions = idx2.query(MARKER);
        std::sort(positions.begin(), positions.end());
        std::vector<uint32_t> expected = {0, 1000};
        assert(positions == expected && "SPEED: query mismatch");
        std::cout << "  query: OK (positions: ";
        for (auto p : positions) std::cout << p << " ";
        std::cout << ")\n";

        // verify_match at each hit: the k-mer itself should produce 0 mismatches.
        const auto packed = pack_read(MARKER);
        for (uint32_t pos : positions) {
            const int mm = idx2.verify_match(pos, packed.data(),
                                             static_cast<uint32_t>(MARKER.size()));
            assert(mm == 0 && "SPEED: verify_match returned non-zero for exact hit");
            std::cout << "  verify_match @ " << pos << ": " << mm << " mismatches — OK\n";
        }

        // Verify a 1-mismatch read: change first base A→C in MARKER.
        std::string one_mm = MARKER;
        one_mm[0] = 'C';  // A→C
        const auto packed_mm = pack_read(one_mm);
        const int mm1 = idx2.verify_match(0, packed_mm.data(),
                                          static_cast<uint32_t>(one_mm.size()));
        assert(mm1 == 1 && "SPEED: verify_match should report 1 mismatch");
        std::cout << "  verify_match 1-mm @ 0: " << mm1 << " mismatch — OK\n";

        // Query wrong-length pattern returns empty.
        assert(idx2.query("ACGT").empty() && "wrong-length query should be empty");
        // Query pattern with N returns empty.
        std::string n_pat = MARKER;
        n_pat[3] = 'N';
        assert(idx2.query(n_pat).empty() && "N-containing query should be empty");
        std::cout << "  edge-case queries: OK\n";
    }

    // -----------------------------------------------------------------------
    // load(MEMORY) + query
    // -----------------------------------------------------------------------
    std::cout << "=== load(MEMORY) ===\n";
    {
        KmerIndex idx3(K);
        idx3.load(idx_path, LoadMode::MEMORY);

        auto positions = idx3.query(MARKER);
        std::sort(positions.begin(), positions.end());
        std::vector<uint32_t> expected = {0, 1000};
        assert(positions == expected && "MEMORY: query mismatch");
        std::cout << "  query: OK (positions: ";
        for (auto p : positions) std::cout << p << " ";
        std::cout << ")\n";

        const auto packed = pack_read(MARKER);
        for (uint32_t pos : positions) {
            const int mm = idx3.verify_match(pos, packed.data(),
                                             static_cast<uint32_t>(MARKER.size()));
            assert(mm == 0 && "MEMORY: verify_match returned non-zero for exact hit");
            std::cout << "  verify_match @ " << pos << ": " << mm << " mismatches — OK\n";
        }
    }

    // -----------------------------------------------------------------------
    // Error handling: bad magic
    // -----------------------------------------------------------------------
    std::cout << "=== error handling ===\n";
    {
        bool threw = false;
        try {
            KmerIndex idx4(K);
            idx4.load(fasta_path, LoadMode::MEMORY);  // FASTA is not an index file
        } catch (const std::runtime_error& e) {
            threw = true;
            std::cout << "  bad-magic exception: " << e.what() << " — OK\n";
        }
        assert(threw && "expected exception on bad magic");
    }

    std::cout << "\nAll tests passed.\n";
    return 0;
}

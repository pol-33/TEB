//   GenomeIndex (abstract base)
//     ├── KmerIndex        — Hash-table of k-mers → position lists
//     ├── SuffixArrayIndex — Suffix array with binary/LCP search
//     └── FMIndex          — BWT-based FM-Index with Occ/C tables
//
// The --optimize-memory flag selects between two internal representations:
//
//   OptimizeMode::SPEED  (default)
//     • Genome stored as flat std::vector<char> or memory-mapped file.
//     • All auxiliary tables (SA, Occ, hash buckets) kept uncompressed in RAM.
//     • Maximises query throughput at the cost of higher RSS.
//
//   OptimizeMode::MEMORY
//     • Genome encoded in 2-bit packing (A=00, C=01, G=10, T=11).
//     • 'N' bases handled via a separate Run-Length Encoded side-table of
//       (start_pos, run_length) pairs, so they cost almost nothing.
//     • Suffix arrays sampled (e.g., every 32nd entry); rank/select structures
//       replace flat Occ tables — trades decode overhead for lower footprint.

#ifndef TEB_INDEX_HPP
#define TEB_INDEX_HPP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

// ---------- enums --------------------------------------------------------- //

/// Controls the internal storage strategy of every index.
enum class OptimizeMode {
    SPEED,   // Flat / uncompressed — fast random access
    MEMORY   // 2-bit packing + RLE for 'N', sampled structures — low RSS
};

/// Indexing algorithm choice, set via --algo on the CLI.
enum class IndexAlgo {
    KMER,    // K-mer hash table
    SA,      // Suffix Array
    FM       // FM-Index (BWT-based)
};

// ---------- abstract base ------------------------------------------------- //

/// Base interface for every genome index backend.
/// Derived classes implement the four core operations and may add
/// algorithm-specific helpers as private methods.
class GenomeIndex {
public:
    virtual ~GenomeIndex() = default;

    /// Parse the FASTA file and construct the in-memory index.
    virtual void build(const std::string& fasta_path) = 0;

    /// Serialise the index to a binary file on disk.
    virtual void save(const std::string& idx_path) const = 0;

    /// Deserialise a previously saved index from disk.
    virtual void load(const std::string& idx_path) = 0;

    /// Return every occurrence position of @p pattern in the indexed genome.
    virtual std::vector<size_t> query(const std::string& pattern) const = 0;

    // ---- factory --------------------------------------------------------- //

    /// Instantiate the right derived class based on algorithm + optimisation.
    static std::unique_ptr<GenomeIndex> create(IndexAlgo algo, OptimizeMode mode);
};

// ---------- K-mer hash index ---------------------------------------------- //

/// Hash-table index: for each k-mer, stores a sorted list of genome positions.
///
/// SPEED mode  — genome kept as std::vector<char>; hash map with open addressing.
/// MEMORY mode — genome in 2-bit packed form; positions delta-encoded + varint.
class KmerIndex : public GenomeIndex {
public:
    explicit KmerIndex(OptimizeMode mode);

    void build(const std::string& fasta_path) override;
    void save(const std::string& idx_path) const override;
    void load(const std::string& idx_path) override;
    std::vector<size_t> query(const std::string& pattern) const override;

private:
    OptimizeMode mode_;
    // TODO(milestone-2): flat genome buffer / 2-bit packed genome
    // TODO(milestone-2): hash table k-mer → position list
};

// ---------- Suffix Array -------------------------------------------------- //

/// Classic suffix array with binary search (or LCP-accelerated search).
///
/// SPEED mode  — full 64-bit SA in RAM; optional LCP array for faster search.
/// MEMORY mode — sampled SA (every s-th entry) + inverse-SA / Ψ function to reconstruct missing entries on the fly.
class SuffixArrayIndex : public GenomeIndex {
public:
    explicit SuffixArrayIndex(OptimizeMode mode);

    void build(const std::string& fasta_path) override;
    void save(const std::string& idx_path) const override;
    void load(const std::string& idx_path) override;
    std::vector<size_t> query(const std::string& pattern) const override;

private:
    OptimizeMode mode_;
    // TODO(milestone-2): SA vector, optional LCP, genome text
};

// ---------- FM-Index ------------------------------------------------------ //

/// BWT-based FM-Index supporting O(m) exact pattern matching via backward search.
///
/// SPEED mode  — flat Occ table (alphabet × n) for O(1) rank queries; full BWT + C array in contiguous memory.
/// MEMORY mode — compressed BWT with run-length encoding; Occ stored as sampled checkpoints (e.g., every 128 rows) + on-the-fly scan.
class FMIndex : public GenomeIndex {
public:
    explicit FMIndex(OptimizeMode mode);

    void build(const std::string& fasta_path) override;
    void save(const std::string& idx_path) const override;
    void load(const std::string& idx_path) override;
    std::vector<size_t> query(const std::string& pattern) const override;

private:
    OptimizeMode mode_;
    // TODO(milestone-2): BWT, C array, Occ table / compressed wavelet tree
};

#endif // TEB_INDEX_HPP

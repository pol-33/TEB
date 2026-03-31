#ifndef MAPPER_MEMORY_FASTA_READER_HPP
#define MAPPER_MEMORY_FASTA_READER_HPP

#include <cstdint>
#include <fstream>
#include <mutex>
#include <string>
#include <vector>

namespace mapper_memory {

struct ChromInfo {
    std::string name;
    uint64_t offset = 0;
    uint64_t length = 0;
};

struct FastaData {
    uint64_t genome_length = 0;
    std::vector<uint8_t> packed_genome;
    std::vector<uint8_t> text;
    std::vector<ChromInfo> chromosomes;

    uint64_t text_length() const {
        return static_cast<uint64_t>(text.size());
    }
};

struct FastaChromosome {
    std::string name;
    uint64_t offset = 0;
    std::string sequence;
};

class FastaReader {
public:
    explicit FastaReader(const std::string& path);

    bool next(FastaChromosome& chrom);

private:
    std::ifstream in_;
    std::string pending_header_;
    uint64_t next_offset_ = 0;
};

FastaData load_fasta(const std::string& path);

// ============================================================================
// Indexed FASTA for streaming reference access (memory-efficient)
// ============================================================================

// FASTA index entry for random access
struct FastaIndexEntry {
    std::string name;
    uint64_t length = 0;           // sequence length in bases
    uint64_t file_offset = 0;      // byte offset of first base in file
    uint32_t line_bases = 0;       // bases per line
    uint32_t line_bytes = 0;       // bytes per line (including newline)
};

// Memory-efficient streaming FASTA reader with indexed random access
// Only loads small regions on demand, does not store full genome
class IndexedFasta {
public:
    IndexedFasta();
    explicit IndexedFasta(const std::string& fasta_path);
    ~IndexedFasta();

    IndexedFasta(const IndexedFasta&) = delete;
    IndexedFasta& operator=(const IndexedFasta&) = delete;
    IndexedFasta(IndexedFasta&& other) noexcept;
    IndexedFasta& operator=(IndexedFasta&& other) noexcept;

    void open(const std::string& fasta_path);
    void close();
    bool is_open() const;

    // Get number of chromosomes/contigs
    std::size_t chromosome_count() const;
    
    // Get chromosome info by index
    const FastaIndexEntry& chromosome(std::size_t index) const;
    
    // Find chromosome by name (returns SIZE_MAX if not found)
    std::size_t find_chromosome(const std::string& name) const;

    // Extract reference sequence at given position
    // Thread-safe: uses mutex internally
    void extract(std::size_t chrom_index, uint64_t pos, uint64_t length, std::string& out);
    
    // Extract reference with chromosome name
    void extract(const std::string& chrom_name, uint64_t pos, uint64_t length, std::string& out);

private:
    void build_index();
    void load_fai_index(const std::string& fai_path);
    void save_fai_index(const std::string& fai_path) const;
    uint64_t file_position_for_base(const FastaIndexEntry& entry, uint64_t base_pos) const;

    std::string fasta_path_;
    std::ifstream file_;
    std::vector<FastaIndexEntry> index_;
    mutable std::mutex mutex_;  // For thread-safe access
};

}  // namespace mapper_memory

#endif

#pragma once

#include "index.hpp"
#include "../genome/genome.hpp"
#include "unordered_dense.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class KmerIndex : public GenomeIndex {
public:
    explicit KmerIndex(size_t k = 14);
    ~KmerIndex();

    void build(const std::string& fasta_path) override;
    void save(const std::string& idx_path) const override;
    void load(const std::string& idx_path,
              LoadMode mode = LoadMode::MEMORY) override;
    std::vector<uint32_t> query(const std::string& pattern) const override;

    int verify_match(uint32_t genome_pos,
                     const uint8_t* read_packed,
                     uint32_t read_len) const;

    uint32_t genome_size() const noexcept { return genome_size_; }
    size_t   k_val()       const noexcept { return k_; }
    std::string chrom_name() const { return chrom_name_; }
    std::string genome_substr(uint32_t pos, uint32_t len) const;
    // Skip k-mers with more than max_freq positions (0 = unlimited).
    void     set_max_freq(uint32_t mf) noexcept { max_freq_ = mf; }

private:
    size_t   k_;
    uint32_t max_freq_ = 0;

    // CSR layout: table_[key] = {start_in_flat, count}.
    // flat_ holds all position lists concatenated.
    ankerl::unordered_dense::map<uint64_t, std::pair<uint32_t, uint32_t>> table_;
    std::vector<uint32_t> flat_;

    // Genome storage (owns data in build mode; empty in load mode).
    GenomeStorage genome_;

    // Raw pointers to packed genome and N-mask.
    // In build mode: point into genome_'s internal storage.
    // In load mode:  point directly into the mmap'd file (zero copy).
    const uint8_t*  packed_ptr_  = nullptr;
    const uint64_t* nmask_ptr_   = nullptr;
    uint32_t        genome_size_ = 0;
    std::string     chrom_name_;

    // mmap tracking (non-zero only after load()).
    void*  mmap_ptr_  = nullptr;
    size_t mmap_size_ = 0;
    int    mmap_fd_   = -1;
};

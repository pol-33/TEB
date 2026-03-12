#pragma once

#include "index.hpp"
#include "unordered_dense.h"
#include <cstddef>
#include <cstdint>
#include <vector>

class KmerIndex : public GenomeIndex {
public:
    explicit KmerIndex(size_t k = 14);

    void build(const std::string& fasta_path) override;
    void save(const std::string& idx_path) const override;
    void load(const std::string& idx_path,
              LoadMode mode = LoadMode::MEMORY) override;
    std::vector<uint32_t> query(const std::string& pattern) const override;

    int verify_match(uint32_t genome_pos,
                     const uint8_t* read_packed,
                     uint32_t read_len) const;

private:
    [[maybe_unused]] size_t k_;
    ankerl::unordered_dense::map<uint64_t, std::vector<uint32_t>> table_;
};

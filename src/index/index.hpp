#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

enum class IndexAlgo    { KMER, SA, FM };
enum class OptimizeMode { SPEED, MEMORY };
enum class LoadMode     { SPEED, MEMORY };

class GenomeIndex {
public:
    virtual void build(const std::string& fasta_path) = 0;
    virtual void save(const std::string& idx_path) const = 0;
    virtual void load(const std::string& idx_path,
                      LoadMode mode = LoadMode::MEMORY) = 0;
    virtual std::vector<uint32_t> query(const std::string& pattern) const = 0;
    virtual ~GenomeIndex() = default;

    static std::unique_ptr<GenomeIndex> create(IndexAlgo, OptimizeMode);
};

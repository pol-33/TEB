#include "index.hpp"

#include <stdexcept>
#include <string>

// ---- factory ------------------------------------------------------------- //

std::unique_ptr<GenomeIndex> GenomeIndex::create(IndexAlgo algo, OptimizeMode mode) {
    switch (algo) {
        case IndexAlgo::KMER: return std::make_unique<KmerIndex>(mode);
        case IndexAlgo::SA:   return std::make_unique<SuffixArrayIndex>(mode);
        case IndexAlgo::FM:   return std::make_unique<FMIndex>(mode);
    }
    throw std::runtime_error("Unknown IndexAlgo value");
}

// ---- KmerIndex stubs ----------------------------------------------------- //

KmerIndex::KmerIndex(OptimizeMode mode) : mode_(mode) {}

void KmerIndex::build(const std::string& /*fasta_path*/) {
    (void)mode_;
    throw std::runtime_error("KmerIndex::build() not yet implemented");
}
void KmerIndex::save(const std::string& /*idx_path*/) const {
    throw std::runtime_error("KmerIndex::save() not yet implemented");
}
void KmerIndex::load(const std::string& /*idx_path*/) {
    throw std::runtime_error("KmerIndex::load() not yet implemented");
}
std::vector<size_t> KmerIndex::query(const std::string& /*pattern*/) const {
    throw std::runtime_error("KmerIndex::query() not yet implemented");
}

// ---- SuffixArrayIndex stubs ---------------------------------------------- //

SuffixArrayIndex::SuffixArrayIndex(OptimizeMode mode) : mode_(mode) {}

void SuffixArrayIndex::build(const std::string& /*fasta_path*/) {
    (void)mode_;
    throw std::runtime_error("SuffixArrayIndex::build() not yet implemented");
}
void SuffixArrayIndex::save(const std::string& /*idx_path*/) const {
    throw std::runtime_error("SuffixArrayIndex::save() not yet implemented");
}
void SuffixArrayIndex::load(const std::string& /*idx_path*/) {
    throw std::runtime_error("SuffixArrayIndex::load() not yet implemented");
}
std::vector<size_t> SuffixArrayIndex::query(const std::string& /*pattern*/) const {
    throw std::runtime_error("SuffixArrayIndex::query() not yet implemented");
}

// ---- FMIndex stubs ------------------------------------------------------- //

FMIndex::FMIndex(OptimizeMode mode) : mode_(mode) {}

void FMIndex::build(const std::string& /*fasta_path*/) {
    (void)mode_;
    throw std::runtime_error("FMIndex::build() not yet implemented");
}
void FMIndex::save(const std::string& /*idx_path*/) const {
    throw std::runtime_error("FMIndex::save() not yet implemented");
}
void FMIndex::load(const std::string& /*idx_path*/) {
    throw std::runtime_error("FMIndex::load() not yet implemented");
}
std::vector<size_t> FMIndex::query(const std::string& /*pattern*/) const {
    throw std::runtime_error("FMIndex::query() not yet implemented");
}

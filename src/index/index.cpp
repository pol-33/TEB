#include "index.hpp"
#include "kmer_index.hpp"
#include <stdexcept>

std::unique_ptr<GenomeIndex> GenomeIndex::create(IndexAlgo algo, OptimizeMode) {
    switch (algo) {
        case IndexAlgo::KMER: return std::make_unique<KmerIndex>();
        default:              throw std::runtime_error("not implemented");
    }
}

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "bwt.hpp"
#include "fasta_reader.hpp"
#include "fm_index.hpp"

namespace {

struct Config {
    std::string reference_path;
    std::string index_path;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-R" && i + 1 < argc) {
            config.reference_path = argv[++i];
        } else if (arg == "-I" && i + 1 < argc) {
            config.index_path = argv[++i];
        } else {
            throw std::runtime_error("usage: indexer -R genome.fa -I genome.idx");
        }
    }
    if (config.reference_path.empty() || config.index_path.empty()) {
        throw std::runtime_error("usage: indexer -R genome.fa -I genome.idx");
    }
    return config;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config config = parse_args(argc, argv);

        std::cerr << "[indexer] loading FASTA: " << config.reference_path << '\n';
        mapper_memory::FastaData fasta = mapper_memory::load_fasta(config.reference_path);
        std::cerr << "[indexer] loaded " << fasta.chromosomes.size() << " chromosome(s), "
                  << fasta.genome.size() << " bases\n";

        std::cerr << "[indexer] building suffix array with SA-IS\n";
        std::vector<uint32_t> suffix_array = mapper_memory::build_suffix_array_sais(fasta.genome);

        std::cerr << "[indexer] deriving packed BWT\n";
        mapper_memory::BWTData bwt = mapper_memory::build_bwt(fasta.genome, suffix_array);
        suffix_array.clear();
        suffix_array.shrink_to_fit();

        std::cerr << "[indexer] building FM-index tables\n";
        mapper_memory::OwnedFMIndex index =
            mapper_memory::build_fm_index(fasta.genome, fasta.chromosomes, bwt, 128, 32);

        bwt.suffix_array.clear();
        bwt.suffix_array.shrink_to_fit();

        std::cerr << "[indexer] writing index: " << config.index_path << '\n';
        const std::size_t bytes_written = index.write(config.index_path);
        std::cerr << "[indexer] final index size: " << bytes_written << " bytes\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}

#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "seed_index.hpp"

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

        std::cerr << "[indexer] building a low-memory " << mapper_memory::kSeedLength
                  << "-mer seed index from: " << config.reference_path << '\n';
        const mapper_memory::IndexBuildStats stats =
            mapper_memory::build_seed_index(config.reference_path, config.index_path);
        std::cerr << "[indexer] indexed " << stats.genome_length << " bases across "
                  << stats.num_chromosomes << " chromosome(s)\n";
        std::cerr << "[indexer] indexed positions: " << stats.indexed_positions << '\n';
        std::cerr << "[indexer] final index size: " << stats.file_size << " bytes\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}

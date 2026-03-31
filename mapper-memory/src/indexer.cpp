#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bwt.hpp"
#include "fasta_reader.hpp"
#include "fm_index.hpp"
#include "memory_stats.hpp"
#include "nucleotide.hpp"

namespace {

struct Config {
    std::string reference_path;
    std::string index_path;
    int opt_level = 0;  // Optimization level 0-5
};

void print_usage() {
    std::cerr << "Usage: indexer -R genome.fa -I genome.idx [-L level]\n"
              << "\n"
              << "Optimization levels:\n"
              << "  0 (baseline): occ=256, sa=32, genome=yes  (largest, fastest)\n"
              << "  1 (no-genome): occ=256, sa=32, genome=no  (no packed genome)\n"
              << "  2 (sparse-sa): occ=256, sa=128, genome=no (sparser SA)\n"
              << "  3 (very-sparse): occ=512, sa=256, genome=no (very sparse)\n"
              << "  4 (ultra-sparse): occ=1024, sa=512, genome=no (ultra sparse)\n"
              << "  5 (extreme): occ=2048, sa=1024, genome=no (extreme sparse, slowest)\n";
}

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-R" && i + 1 < argc) {
            config.reference_path = argv[++i];
        } else if (arg == "-I" && i + 1 < argc) {
            config.index_path = argv[++i];
        } else if (arg == "-L" && i + 1 < argc) {
            config.opt_level = std::stoi(argv[++i]);
            if (config.opt_level < 0 || config.opt_level > 5) {
                throw std::runtime_error("optimization level must be 0-5");
            }
        } else if (arg == "-h" || arg == "--help") {
            print_usage();
            std::exit(0);
        } else {
            print_usage();
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (config.reference_path.empty() || config.index_path.empty()) {
        print_usage();
        throw std::runtime_error("missing required arguments");
    }
    return config;
}

mapper_memory::IndexConfig get_index_config(int level) {
    switch (level) {
        case 0: return mapper_memory::IndexConfig::baseline();
        case 1: return mapper_memory::IndexConfig::level1_no_genome();
        case 2: return mapper_memory::IndexConfig::level2_sparse_sa();
        case 3: return mapper_memory::IndexConfig::level3_very_sparse();
        case 4: return mapper_memory::IndexConfig::level4_ultra_sparse();
        case 5: return mapper_memory::IndexConfig::level5_extreme();
        default: return mapper_memory::IndexConfig::baseline();
    }
}

std::vector<uint8_t> sequence_to_text(const std::string& sequence) {
    std::vector<uint8_t> text;
    text.reserve(sequence.size() + 1U);
    for (char base : sequence) {
        text.push_back(mapper_memory::rank_from_base(base));
    }
    text.push_back(mapper_memory::kSentinelRank);
    return text;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config config = parse_args(argc, argv);
        const mapper_memory::IndexConfig idx_config = get_index_config(config.opt_level);

        std::cerr << "[indexer] optimization level " << config.opt_level << ": "
                  << "occ_sample=" << idx_config.occ_sample
                  << ", sa_sample=" << idx_config.sa_sample
                  << ", store_genome=" << (idx_config.store_genome ? "yes" : "no")
                  << "\n";

        mapper_memory::FastaReader reader(config.reference_path);
        mapper_memory::FastaChromosome chrom;
        std::vector<mapper_memory::OwnedChromosomeIndex> chromosomes;
        chromosomes.reserve(256);

        uint64_t total_bases = 0;
        while (reader.next(chrom)) {
            total_bases += chrom.sequence.size();
            std::cerr << "Building chromosome " << chrom.name << " (" << chrom.sequence.size() << " bp)\n";

            std::vector<uint8_t> text = sequence_to_text(chrom.sequence);
            std::vector<uint32_t> suffix_array = mapper_memory::build_suffix_array_sais(text);
            const mapper_memory::BWTData bwt = mapper_memory::build_bwt(text, suffix_array);
            chromosomes.push_back(
                mapper_memory::build_chromosome_index(chrom.name, chrom.sequence, bwt, suffix_array, idx_config));

            suffix_array.clear();
            suffix_array.shrink_to_fit();
            text.clear();
            text.shrink_to_fit();
            chrom.sequence.clear();
            chrom.sequence.shrink_to_fit();

            const mapper_memory::MemoryStats stats = mapper_memory::read_memory_stats();
            std::cerr << "Finished " << chrom.name
                      << " - current RSS " << mapper_memory::bytes_to_mebibytes(stats.current_rss_bytes)
                      << " MiB, peak RSS " << mapper_memory::bytes_to_mebibytes(stats.peak_rss_bytes) << " MiB\n";
        }

        std::cerr << "Genome loaded: " << total_bases << " bases across " << chromosomes.size() << " chromosomes\n";
        const std::size_t file_size = mapper_memory::FMIndexView::write(config.index_path, chromosomes, idx_config);
        std::cerr << "Index written: "
                  << static_cast<double>(file_size) / (1024.0 * 1024.0 * 1024.0)
                  << " GiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}

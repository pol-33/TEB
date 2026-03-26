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
                mapper_memory::build_chromosome_index(chrom.name, chrom.sequence, bwt, suffix_array, 256U, 32U));

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
        const std::size_t file_size = mapper_memory::FMIndexView::write(config.index_path, chromosomes, 256U, 32U);
        std::cerr << "Index written: "
                  << static_cast<double>(file_size) / (1024.0 * 1024.0 * 1024.0)
                  << " GiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}

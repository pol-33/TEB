// Usage:
//   indexer -R <genome.fa> -I <genome.idx> [--algorithm kmer|sa|fm] [--optimize-memory]
//
// Builds a genome index from a FASTA reference and writes it to disk.
// The --algorithm flag selects the indexing backend algorithm (default: fm).
// The --optimize-memory flag switches from speed-optimised (default) to
// memory-optimised internal representations (2-bit packing, sampled SA, etc.).

#include "index.hpp"
#include "kmer_index.hpp"

#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <iostream>
#include <string>

// ---- CLI config ---------------------------------------------------------- //

struct IndexerConfig {
    std::string fasta_path;              // -R <genome.fa>
    std::string idx_path;                // -I <genome.idx>
    IndexAlgo   algorithm = IndexAlgo::FM;   // --algorithm (default: fm)
    OptimizeMode mode = OptimizeMode::SPEED; // --optimize-memory flips to MEMORY
    uint32_t    max_freq  = 0;               // --max-freq N (0 = disabled)
};

// ---- usage / help -------------------------------------------------------- //

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " -R <genome.fa> -I <genome.idx> [options]\n"
        << "\n"
        << "Required arguments:\n"
        << "  -R <path>          Input FASTA reference genome\n"
        << "  -I <path>          Output index file path\n"
        << "\n"
        << "Options:\n"
        << "  --algorithm <kmer|sa|fm>   Indexing algorithm (default: fm)\n"
        << "  --optimize-memory         Optimise for memory instead of speed\n"
        << "  --max-freq <N>            Skip k-mers with >N occurrences (kmer only)\n"
        << "  -h, --help                Show this help message\n";
}

// ---- argument parsing (getopt_long) -------------------------------------- //

// Long-option IDs that don't clash with short-option chars.
enum LongOpt : int {
    OPT_ALGORITHM = 256,
    OPT_OPTMEM    = 257,
    OPT_MAXFREQ   = 258
};

static IndexerConfig parse_args(int argc, char* argv[]) {
    IndexerConfig cfg;

    // Short options: -R, -I, -h
    static const char* short_opts = "R:I:h";

    // Long options table
    static struct option long_opts[] = {
        {"algorithm",       required_argument, nullptr, OPT_ALGORITHM},
        {"optimize-memory", no_argument,       nullptr, OPT_OPTMEM},
        {"max-freq",        required_argument, nullptr, OPT_MAXFREQ},
        {"help",            no_argument,       nullptr, 'h'},
        {nullptr,           0,                 nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'R':
                cfg.fasta_path = optarg;
                break;
            case 'I':
                cfg.idx_path = optarg;
                break;
            case OPT_ALGORITHM:
                if (std::strcmp(optarg, "kmer") == 0)      cfg.algorithm = IndexAlgo::KMER;
                else if (std::strcmp(optarg, "sa") == 0)    cfg.algorithm = IndexAlgo::SA;
                else if (std::strcmp(optarg, "fm") == 0)    cfg.algorithm = IndexAlgo::FM;
                else {
                    std::cerr << "Error: unknown algorithm '" << optarg
                              << "'. Choose from: kmer, sa, fm.\n";
                    std::exit(EXIT_FAILURE);
                }
                break;
            case OPT_OPTMEM:
                cfg.mode = OptimizeMode::MEMORY;
                break;
            case OPT_MAXFREQ: {
                const long v = std::strtol(optarg, nullptr, 10);
                if (v <= 0) {
                    std::cerr << "Error: --max-freq must be a positive integer.\n";
                    std::exit(EXIT_FAILURE);
                }
                cfg.max_freq = static_cast<uint32_t>(v);
                break;
            }
            case 'h':
                print_usage(argv[0]);
                std::exit(EXIT_SUCCESS);
            default:
                print_usage(argv[0]);
                std::exit(EXIT_FAILURE);
        }
    }

    // Validate required arguments
    if (cfg.fasta_path.empty() || cfg.idx_path.empty()) {
        std::cerr << "Error: -R and -I are required.\n\n";
        print_usage(argv[0]);
        std::exit(EXIT_FAILURE);
    }

    return cfg;
}

// ---- helpers ------------------------------------------------------------- //

static const char* algo_name(IndexAlgo a) {
    switch (a) {
        case IndexAlgo::KMER: return "kmer";
        case IndexAlgo::SA:   return "sa";
        case IndexAlgo::FM:   return "fm";
    }
    return "unknown";
}

static const char* mode_name(OptimizeMode m) {
    return m == OptimizeMode::MEMORY ? "memory" : "speed";
}

// ---- main ---------------------------------------------------------------- //

int main(int argc, char* argv[]) {
    IndexerConfig cfg = parse_args(argc, argv);

    std::cerr << "[indexer] Reference : " << cfg.fasta_path << "\n"
              << "[indexer] Output    : " << cfg.idx_path    << "\n"
              << "[indexer] Algorithm : " << algo_name(cfg.algorithm) << "\n"
              << "[indexer] Optimise  : " << mode_name(cfg.mode) << "\n";

    // Create the appropriate index backend via the index factory.
    auto index = GenomeIndex::create(cfg.algorithm, cfg.mode);

    // Apply max-freq filter before building (kmer index only).
    if (cfg.max_freq > 0) {
        if (auto* ki = dynamic_cast<KmerIndex*>(index.get()))
            ki->set_max_freq(cfg.max_freq);
        else
            std::cerr << "[indexer] Warning: --max-freq ignored for non-kmer algorithm.\n";
    }

    // Build the index from the FASTA reference.
    std::cerr << "[indexer] Building index...\n";
    index->build(cfg.fasta_path);

    // Write the index to disk.
    std::cerr << "[indexer] Saving index to " << cfg.idx_path << "...\n";
    index->save(cfg.idx_path);

    std::cerr << "[indexer] Done.\n";
    return EXIT_SUCCESS;
}

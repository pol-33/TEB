// Usage:
//   mapper -I <genome.idx> -i <reads.fastq> -o <output.sam> -k <max-errors>
//
// Loads a pre-built genome index, maps each read from the FASTQ file allowing
// up to k mismatches/edits, and writes alignments in SAM format.

#include "index.hpp"

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>

// ---- CLI config ---------------------------------------------------------- //

struct MapperConfig {
    std::string idx_path;       // -I <genome.idx>
    std::string reads_path;     // -i <reads.fastq>
    std::string output_path;    // -o <output.sam>
    int         max_errors = 0; // -k <max-errors>
};

// ---- usage / help -------------------------------------------------------- //

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -I <genome.idx> -i <reads.fastq> -o <output.sam> -k <max-errors>\n"
        << "\n"
        << "Required arguments:\n"
        << "  -I <path>   Pre-built genome index file\n"
        << "  -i <path>   Input FASTQ reads file\n"
        << "  -o <path>   Output SAM alignment file\n"
        << "  -k <int>    Maximum number of errors (mismatches/edits)\n"
        << "\n"
        << "Options:\n"
        << "  -h, --help  Show this help message\n";
}

// ---- argument parsing (getopt_long) -------------------------------------- //

static MapperConfig parse_args(int argc, char* argv[]) {
    MapperConfig cfg;

    static const char* short_opts = "I:i:o:k:h";

    static struct option long_opts[] = {
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'I':
                cfg.idx_path = optarg;
                break;
            case 'i':
                cfg.reads_path = optarg;
                break;
            case 'o':
                cfg.output_path = optarg;
                break;
            case 'k': {
                char* endptr = nullptr;
                long val = std::strtol(optarg, &endptr, 10);
                if (endptr == optarg || *endptr != '\0' || val < 0) {
                    std::cerr << "Error: -k requires a non-negative integer.\n";
                    std::exit(EXIT_FAILURE);
                }
                cfg.max_errors = static_cast<int>(val);
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

    // Validate all four required arguments
    if (cfg.idx_path.empty() || cfg.reads_path.empty() ||
        cfg.output_path.empty()) {
        std::cerr << "Error: -I, -i, -o, and -k are all required.\n\n";
        print_usage(argv[0]);
        std::exit(EXIT_FAILURE);
    }

    return cfg;
}

// ---- main ---------------------------------------------------------------- //

int main(int argc, char* argv[]) {
    MapperConfig cfg = parse_args(argc, argv);

    std::cerr << "[mapper] Index      : " << cfg.idx_path    << "\n"
              << "[mapper] Reads      : " << cfg.reads_path  << "\n"
              << "[mapper] Output     : " << cfg.output_path << "\n"
              << "[mapper] Max errors : " << cfg.max_errors  << "\n";

    // TODO(milestone-3): detect which algorithm was used from the index header,
    //   load the index via GenomeIndex::load(), read FASTQ, align, write SAM.
    std::cerr << "[mapper] Mapping logic not yet implemented.\n";

    return EXIT_SUCCESS;
}

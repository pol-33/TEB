#include "fasta_parser.hpp"
#include "fastq_parser.hpp"

#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

using namespace std;

enum class Format {
    FASTA,
    FASTQ
};

struct Config {
    string input;
    string output;
    Format format;
    int qmin = 0;
    unsigned int kmer_length = 0;
    bool per_seq_stats = false;
    bool low_mem = false;
    size_t stream_buf_bytes = 3UL * 1024UL * 1024UL;  // default 3 MB
};

void usage() {
    cout << "./teb <args> [options]       \n";
    cout << "\t Compulsory args:          \n";
    cout << "\t\t -i <filename>     Input filename  \n";
    cout << "\t\t -f <fasta/fastq>  Format of the input filename\n";
    cout << "\t Optional args:            \n";
    cout << "\t\t -o <filename>     Output filename (omit for stats-only mode) \n";
    cout << "\t\t -qmin <int value> Minimum Quality for bases (FASTQ only)\n";
    cout << "\t\t -k <int value>    Kmer length for preprocessing text\n";
    cout << "\t\t -s                Print per-sequence statistics (off by default)\n";
    cout << "\t\t -m [MB]           Low-memory streaming mode (default 3 MB buffer)\n";
    cout << "\t\t                   e.g. -m 0.5, -m 3, -m 64  (float MB accepted)\n";
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    bool format_set = false;

    for (int i = 1; i < argc; i++) {
        string arg = argv[i];

        if (arg == "-i") {
            cfg.input = argv[++i];
        } else if (arg == "-o") {
            cfg.output = argv[++i];
        } else if (arg == "-f") {
            string f = argv[++i];

            if (f == "fasta") cfg.format = Format::FASTA;
            else if (f == "fastq") cfg.format = Format::FASTQ;
            else throw runtime_error("[parse_args] Invalid format");

            format_set = true;
        } else if (arg == "-qmin") {
            cfg.qmin = stoi(argv[++i]);
        } else if (arg == "-k") {
            cfg.kmer_length = stoi(argv[++i]);
        } else if (arg == "-s") {
            cfg.per_seq_stats = true;
        } else if (arg == "-m") {
            cfg.low_mem = true;
            // Optional next arg: buffer size in MB (float accepted)
            if (i + 1 < argc) {
                char* endptr;
                double mb = strtod(argv[i + 1], &endptr);
                if (endptr != argv[i + 1] && *endptr == '\0' && mb > 0.0) {
                    cfg.stream_buf_bytes = (size_t)(mb * 1024.0 * 1024.0);
                    ++i;
                }
            }
        } else throw runtime_error("[parse_args] Unknown argument: " + arg);
    }

    if (cfg.input.empty() || !format_set)
        throw runtime_error("[parse_args] Missing required arguments");

    return cfg;
}

int main(int argc, char* argv[]) {
    clock_t t = clock();
    try {
        Config cfg = parse_args(argc, argv);

        if (cfg.format == Format::FASTA) fasta_parser(cfg.input, cfg.output, cfg.kmer_length, cfg.per_seq_stats, cfg.low_mem, cfg.stream_buf_bytes);
        else if (cfg.format == Format::FASTQ) fastq_parser(cfg.input, cfg.output, cfg.qmin, cfg.per_seq_stats, cfg.low_mem, cfg.stream_buf_bytes);

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        usage();
        return 1;
    }

    t = clock() -t;
    long seconds = 0;
    long remainder = t;

    /* Count seconds without division */
    while (remainder >= 1000000) {
        remainder -= 1000000;
        seconds++;
    }

    /* Print as seconds.microseconds */
    printf("Total time std (seconds): %ld.%06ld\n", seconds, remainder);
    return 0;
}
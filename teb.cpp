#include "fasta_parser.hpp"

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
};

void usage() {
    cout << "./teb <args> [options]       \n";
    cout << "\t Compulsory args:          \n";
    cout << "\t\t -i <filename>    Input filename  \n";
    cout << "\t\t -o <filename>    Output filename \n";
    cout << "\t\t -f <fasta/fastq> Format of the input filename \n";
    cout << "\t Optional args:            \n";
    cout << "\t\t -qmin <int value>      \n";
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
        } else throw runtime_error("[parse_args] Unknown argument: " + arg);
    }

    if (cfg.input.empty() || cfg.output.empty() || !format_set)
        throw runtime_error("[parse_args] Missing required arguments");

    return cfg;
}

int main(int argc, char* argv[]) {
    clock_t t = clock();
    try {
        Config cfg = parse_args(argc, argv);

        if (cfg.format == Format::FASTA) fasta_parser(cfg.input, cfg.output);
        else if (cfg.format == Format::FASTA) fasta_parser(cfg.input, cfg.output);

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        usage();
        return 1;
    }

    t = clock() -t;
    printf("Total time std (micro-seconds): %ld \n", t);
    return 0;
}
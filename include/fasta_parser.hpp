#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "common.hpp"

#include "utils.hpp"

using namespace std;

// Estructura per guardar cada registre FASTA
struct FastaRecord {
    string header;
    string sequence;
};

int fasta_parser(const string& input_file, const string& output_file = "", const unsigned int kmer_length = 0, bool per_seq_stats = false, bool low_mem = true, size_t stream_buf = (size_t)(0.5 * 1024UL * 1024UL));

#endif
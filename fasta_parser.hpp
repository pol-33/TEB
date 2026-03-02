#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>

#include "utils.hpp"

using namespace std;

// Estructura per guardar cada registre FASTA
struct FastaRecord {
    string header;
    string sequence;
};

struct SequenceStats {
    long long length;
    long long gc_count;
    double gc_content;
};

struct GlobalStats {
    int num_sequences;
    long long total_length;
    long long total_gc_count;
    double overall_gc_content;
};

int fasta_parser(const string& input_file, const string& output_file, const unsigned int kmer_length);

#endif
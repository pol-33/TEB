#ifndef __FASTQ_PARSER__
#define __FASTQ_PARSER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "utils.hpp"

using namespace std;

// Estructura per guardar cada registre FASTA
struct FastQRecord {
    string header;
    string quality_sequence;
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
    map<string, SequenceStats> per_sequence;
    long long minimum;
    long long maximum;
    long long avg_read_len;
};

int fasta_parser(const string& input_file, const string& output_file, const int qmin);

#endif
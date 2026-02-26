#ifndef __FASTQ_PARSER__
#define __FASTQ_PARSER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "common.hpp"

using namespace std;

// Estructura per guardar cada registre FASTQ
struct FastqRecord {
    string header;
    string sequence;
    string quality_sequence;
};

struct FastqGlobalStats : public GlobalStats {
    long long minimum;
    long long maximum;
    long long avg_read_len;
};

int fastq_parser(const string& input_file, const string& output_file, const int qmin);

#endif

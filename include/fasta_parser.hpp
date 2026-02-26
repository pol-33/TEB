#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include "common.hpp"

using namespace std;

// Estructura per guardar cada registre FASTA
struct FastaRecord {
    string header;
    string sequence;
};

int fasta_parser(const string& input_file, const string& output_file = "");

#endif

#include "fasta_parser.hpp"

#define KMER_INDEX_SIZE 10

using namespace std;

// Compute statistics
static inline void computeStatistics(GlobalStats& gStats, const string& sequence) {
    int chars_read = sequence.length();
    gStats.num_sequences += 1;
    gStats.total_length += chars_read;

    int local_gc_count = 0;
    for (char c : sequence) {
        local_gc_count += gc_matching(c);
    }

    gStats.total_gc_count += local_gc_count;

    // Calculate global GC
    gStats.overall_gc_content = (chars_read > 0)
                        ? (double)local_gc_count / chars_read * 100.0
                        : 0.0;
    return;
}

// Parsing of the FASTA file
void parseFastaFile(const string& infile, const string& outfile, GlobalStats& stats, const int kmer_length) {
    ifstream in(infile);
    ofstream out(outfile);
    kmer_table_t kmer_indxs;

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    } else if (!out.is_open()) {
        cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
        return;
    }

    string line;
    string current_header = "";

    while (getline(in, line)) {
        if (line.empty()) continue;

        // Character '>', indicates new header
        if (line[0] == '>') {
            current_header = line.substr(1);
            out << current_header << "\n";
        } else {
            // Otherwise: Sequence line -> concat
            out << line;
            computeStatistics(stats, line);
            update_kmer_table(line, kmer_indxs, kmer_length);
        }
    }

    #ifdef DEBUG
    print_kmer_table(kmer_indxs);
    #endif

    in.close();
    out.close();
    return;
}

// Print the statistics
void printStatistics(const GlobalStats& stats) {
    cout << "============================================" << endl;
    cout << "       SUMMARY GENOMIC ANALYSIS" << endl;
    cout << "============================================" << endl;
    cout << "GLOBAL SUMMARY:" << endl;
    cout << "  > Total length: " << stats.total_length << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Number of processed sequences: " << stats.num_sequences << endl;
    cout << "============================================" << endl;
}

int fasta_parser(const string& input_file, const string& output_file, const int kmer_length) {

    GlobalStats stats;
    stats.num_sequences = 0;
    stats.total_length = 0;
    stats.total_gc_count = 0;

    // Step 1: READING
    parseFastaFile(input_file, output_file, stats, kmer_length);
    if (stats.num_sequences == 0) {
        cout << "No sequences found." << endl;
        return 0;
    }

    // Step 2: PRINTING STATS
    printStatistics(stats);

    return 0;
}

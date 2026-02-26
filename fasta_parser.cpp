#include "fasta_parser.hpp"

using namespace std;

// Compute statistics
static inline void computeStatistics(GlobalStats& gStats, const string& sequence) {
    int chars_read = sequence.length();
    gStats.num_sequences += 1;
    gStats.total_length += chars_read;

    int local_gc_count = 0;
    for (char c : sequence) {
        local_gc_count+=gc_matching(c);
    }

    gStats.total_gc_count += local_gc_count;

    // Calculate global GC
    gStats.overall_gc_content = (chars_read > 0)
                        ? (double)local_gc_count / chars_read * 100.0
                        : 0.0;
    return;
}

// Parsing of the FASTA file
void parseFastaFile(const string& infile, const string& outfile, GlobalStats& stats) {
    ifstream in(infile);
    ofstream out(outfile);

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    } else if (!out.is_open()) {
        cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
        return;
    }

    string line;
    string current_header = "";
    string current_sequence = "";

    while (getline(in, line)) {
        if (line.empty()) continue;

        // Character '>', indicates new header
        if (line[0] == '>') {
            current_header = line.substr(1);
            out << current_header << "\n";
        } else {
            // Otherwise: Sequence line -> concat
            current_sequence += line;
            out << current_sequence;
            computeStatistics(stats, current_sequence);
        }
    }

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
    cout << "  > Total lenght: " << stats.total_length << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Number of processed sequences: " << stats.num_sequences << endl;
    cout << "============================================" << endl;
}

int fasta_parser(const string& input_file, const string& output_file) {

    GlobalStats stats;
    stats.num_sequences = 0;
    stats.total_length = 0;
    stats.total_gc_count = 0;

    // Step 1: READING
    parseFastaFile(input_file, output_file, stats);
    if (stats.num_sequences == 0) {
        cout << "No sequences found." << endl;
        return 0;
    }

    // Step 2: PRINTING STATS
    printStatistics(stats);

    return 0;
}

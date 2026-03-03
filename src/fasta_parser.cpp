#include "fasta_parser.hpp"

#define KMER_INDEX_SIZE 10

using namespace std;

// Compute statistics
static inline void computeStatistics(GlobalStats& gStats, SequenceStats& seqStats, const string& sequence) {
    long long chars_read = sequence.length();
    gStats.total_length += chars_read;
    seqStats.length    += chars_read;

    int local_gc_count = 0;
    for (char c : sequence) {
        local_gc_count += gc_matching(c);
    }

    gStats.total_gc_count += local_gc_count;
    seqStats.gc_count     += local_gc_count;

    // Calculate global GC from cumulative totals (not just this line)
    gStats.overall_gc_content = (gStats.total_length > 0)
                        ? (double)gStats.total_gc_count / gStats.total_length * 100.0
                        : 0.0;
    return;
}

// Parsing of the FASTA file
void parseFastaFile(const string& infile, const string& outfile, GlobalStats& stats, const unsigned int kmer_length) {
    ifstream in(infile);

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    }

    ofstream out;
    kmer_table_t kmer_indxs;
    
    if (!outfile.empty()) {
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
            return;
        }
    }

    string line;
    string current_header = "";
    SequenceStats cur_seq = {0, 0, 0.0};

    // Helper: finalise current sequence into the per_sequence map
    auto finalizeSequence = [&]() {
        if (!current_header.empty()) {
            cur_seq.gc_content = (cur_seq.length > 0)
                ? (double)cur_seq.gc_count / cur_seq.length * 100.0
                : 0.0;
            stats.per_sequence[current_header] = cur_seq;
        }
    };

    while (getline(in, line)) {
        if (line.empty()) continue;

        // Character '>', indicates new header
        if (line[0] == '>') {
            finalizeSequence();
            current_header = line.substr(1);
            cur_seq = {0, 0, 0.0};
            stats.num_sequences += 1;
            out << current_header << "\n";
        } else {
            // Otherwise: Sequence line -> concat
            out << line;
            computeStatistics(stats, cur_seq, line);
            if (kmer_length > 0) update_kmer_table(line, kmer_indxs, kmer_length);
        }
    }
    finalizeSequence(); // store the last sequence

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
    cout << "Number of processed sequences: " << stats.num_sequences << endl;
    cout << "--------------------------------------------" << endl;
    for (const auto& [id, seq] : stats.per_sequence) {
        cout << "ID: " << id << endl;
        cout << "  > Length: " << seq.length << " bp" << endl;
        cout << "  > GC content: " << fixed << setprecision(2) << seq.gc_content << "%" << endl;
        cout << endl;
        cout << "--------------------------------------------" << endl;
    }
    cout << "GLOBAL SUMMARY:" << endl;
    cout << "  > Total length: " << stats.total_length << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "============================================" << endl;
}

int fasta_parser(const string& input_file, const string& output_file, const unsigned int kmer_length) {

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

#include "fasta_parser.hpp"

#define KMER_INDEX_SIZE 10

using namespace std;

// Compute statistics
static inline void computeStatistics(GlobalStats& gStats, SequenceStats& seqStats, const string& sequence) {
    long long chars_read = sequence.length();
    gStats.total_length += chars_read;
    seqStats.length    += chars_read;

    long long local_gc = count_gc_bulk(sequence.data(), (size_t)chars_read);
    gStats.total_gc_count += local_gc;
    seqStats.gc_count     += local_gc;

    // Calculate global GC from cumulative totals (not just this line)
    gStats.overall_gc_content = (gStats.total_length > 0)
                        ? (double)gStats.total_gc_count / gStats.total_length * 100.0
                        : 0.0;
    return;
}

// Parsing of the FASTA file
void parseFastaFile(const string& infile, const string& outfile, GlobalStats& stats, const unsigned int kmer_length) {
    // Set large buffer BEFORE open() so the streambuf uses it from the first read
    static char in_buf[IO_BUFFER_SIZE];
    ifstream in;
    in.rdbuf()->pubsetbuf(in_buf, sizeof(in_buf));
    in.open(infile);

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    }

    static char out_buf[IO_BUFFER_SIZE];
    ofstream out;
    kmer_table_t kmer_indxs;
    
    if (!outfile.empty()) {
        out.rdbuf()->pubsetbuf(out_buf, sizeof(out_buf));
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
            if (out.is_open()) {
                out.put('>');
                out.write(current_header.data(), (streamsize)current_header.size());
                out.put('\n');
            }
        } else {
            // Sequence line: raw write avoids format overhead
            if (out.is_open()) {
                out.write(line.data(), (streamsize)line.size());
                out.put('\n');
            }
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

#include "fastq_parser.hpp"
#include "utils.hpp"

#include <climits>
#include <iomanip>

using namespace std;

// Single-pass streaming parser that reads one record at a time and computes stats on-the-fly
static void parseFastqFile(const string& infile, const string& outfile, FastqGlobalStats& stats, const int qmin) {
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
    if (!outfile.empty()) {
        out.rdbuf()->pubsetbuf(out_buf, sizeof(out_buf));
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
            return;
        }
    }

    stats.num_sequences  = 0;
    stats.total_length   = 0;
    stats.total_gc_count = 0;
    stats.minimum        = LLONG_MAX;
    stats.maximum        = 0LL;

    string header, sequence, separator, quality;

    while (getline(in, header)) {
        if (header.empty()) continue;
        if (header[0] != '@') continue;

        if (!getline(in, sequence))  break;
        if (!getline(in, separator)) break;
        if (!getline(in, quality))   break;

        if (sequence.size() != quality.size())
            throw runtime_error("[parseFastqFile] sequence length != quality sequence length");

        // In-place quality trimming using SWAR vectorised scan (see utils.hpp)
        if (qmin > 0) {
            size_t trim = get_trim_limit(quality, qmin);
            sequence.resize(trim);
            quality.resize(trim);
        }

        // Write full FASTQ record with raw write() to avoid << format overhead
        if (out.is_open()) {
            out.write(header.data(),    (streamsize)header.size());    out.put('\n');
            out.write(sequence.data(),  (streamsize)sequence.size());  out.put('\n');
            out.write(separator.data(), (streamsize)separator.size()); out.put('\n');
            out.write(quality.data(),   (streamsize)quality.size());   out.put('\n');
        }

        // Bulk GC count
        long long len = (long long)sequence.size();
        long long gc  = count_gc_bulk(sequence.data(), (size_t)len);

        SequenceStats sStats;
        sStats.length     = len;
        sStats.gc_count   = gc;
        sStats.gc_content = (len > 0) ? (double)gc / len * 100.0 : 0.0;

        stats.per_sequence.emplace(header.substr(1), sStats);
        stats.total_length   += len;
        stats.total_gc_count += gc;
        if (len < stats.minimum) stats.minimum = len;
        if (len > stats.maximum) stats.maximum = len;
        ++stats.num_sequences;
    }

    stats.overall_gc_content = (stats.total_length > 0)
        ? (double)stats.total_gc_count / stats.total_length * 100.0
        : 0.0;
    stats.avg_read_len = (stats.num_sequences > 0)
        ? stats.total_length / stats.num_sequences : 0;
}

// Print the statistics
void printStatistics(const FastqGlobalStats& stats) {
    cout << "============================================" << endl;
    cout << "       SUMMARY GENOMIC ANALYSIS" << endl;
    cout << "============================================" << endl;
    cout << "Number of processed sequences: " << stats.num_sequences << endl;
    cout << "--------------------------------------------" << endl;

    // Iterem sobre el map per mostrar cada seqüència
    // it.first és la clau (header), it.second és l'objecte SequenceStats
    for (auto const& [header, sStats] : stats.per_sequence) {
        cout << "ID: " << header << endl;
        cout << "  > Length: " << sStats.length << " bp" << endl;
        cout << "  > GC content: " << fixed << setprecision(2) << sStats.gc_content << "%" << endl;
        cout << endl;
    }

    cout << "--------------------------------------------" << endl;
    cout << "GLOBAL SUMMARY:" << endl;
    cout << "  > Total lenght: " << stats.total_length << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Minimum read: " << stats.minimum << " bp" << endl;
    cout << "  > Maximum read: " << stats.maximum << " bp" << endl;
    cout << "  > Avg length read: " << stats.avg_read_len << " bp" << endl;
    cout << "============================================" << endl;
}

int fastq_parser(const string& input_file, const string& output_file, const int qmin) {

    FastqGlobalStats stats;

    // Single-pass to parse and compute stats simultaneously
    parseFastqFile(input_file, output_file, stats, qmin);

    if (stats.num_sequences == 0) {
        cout << "No sequences found." << endl;
        return 0;
    }

    printStatistics(stats);
    return 0;
}

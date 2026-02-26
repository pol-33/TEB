#include "fastq_parser.hpp"

using namespace std;

// Parsing of the FASTQ file
vector<FastqRecord> parseFastqFile(const string& infile, const string& outfile, const int qmin = 0) {
    vector<FastqRecord> records;
    ifstream in(infile);

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return records;
    }

    ofstream out;
    if (!outfile.empty()) {
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
            return records;
        }
    }

    string line;
    string current_header = "";
    string current_sequence = "";
    string sequence_quality = "";

    while (getline(in, line)) {
        if (line.empty()) continue;

        // Character '@', indicates new header
        if (line[0] == '@') {
            // save previous content
            if (!current_header.empty()) {
                FastqRecord rec;
                rec.header = current_header;
                rec.sequence = current_sequence;
                rec.quality_sequence = sequence_quality;
                records.push_back(rec);
            }

            current_header = line.substr(1);
            current_sequence = "";
            sequence_quality = "";
            out << current_header << "\n";
        } else {
            current_sequence = line;
            getline(in, line); // "+" separator line
            getline(in, line); // quality line
            sequence_quality = line;
            if (current_sequence.length() != sequence_quality.length()) {
                throw runtime_error("[parseFastqFile] sequence length != quality sequence length");
            }
            if (qmin > 0) {
                // Quality trimming from the right
                int last_good_pos = -1;
                for (int i = (int)current_sequence.length() - 1; i >= 0; i--) {
                    int score = (unsigned char)sequence_quality[i] - 33;
                    if (score >= qmin) {
                        last_good_pos = i;
                        break;
                    }
                }
                if (last_good_pos >= 0) {
                    current_sequence = current_sequence.substr(0, last_good_pos + 1);
                    sequence_quality = sequence_quality.substr(0, last_good_pos + 1);
                } else {
                    current_sequence = "";
                    sequence_quality = "";
                }
            }
        }
    }

    // We save the last record if it exists after the loop ends
    if (!current_header.empty()) {
        FastqRecord rec;
        rec.header = current_header;
        rec.sequence = current_sequence;
        rec.quality_sequence = sequence_quality;
        records.push_back(rec);
    }

    in.close();
    out.close();
    return records;
}

// Compute statistics
FastqGlobalStats computeStatistics(const vector<FastqRecord>& records) {
    FastqGlobalStats gStats;
    gStats.num_sequences = records.size();
    gStats.total_length = 0;
    gStats.total_gc_count = 0;
    gStats.minimum = LLONG_MAX;
    gStats.maximum = 0;

    // Calculate statistics for every registered sequence
    for (const auto& rec : records) {
        SequenceStats sStats;
        sStats.length = rec.sequence.length();
        sStats.gc_count = 0;

        for (char c : rec.sequence) {
            char base = toupper(c);
            if (base == 'G' || base == 'C') {
                sStats.gc_count++;
            }
        }

        // Calculate individual GC
        sStats.gc_content = (sStats.length > 0)
                            ? (double)sStats.gc_count / sStats.length * 100.0
                            : 0.0;

        // Save map and update stats
        gStats.per_sequence[rec.header] = sStats;
        gStats.total_length += sStats.length;
        gStats.total_gc_count += sStats.gc_count;

        // Track min and max read lengths
        if (sStats.length < gStats.minimum) gStats.minimum = sStats.length;
        if (sStats.length > gStats.maximum) gStats.maximum = sStats.length;
    }

    // Calculate global GC
    gStats.overall_gc_content = (gStats.total_length > 0)
                                ? (double)gStats.total_gc_count / gStats.total_length * 100.0
                                : 0.0;

    // Compute average read length
    gStats.avg_read_len = (gStats.num_sequences > 0)
                          ? gStats.total_length / gStats.num_sequences
                          : 0;

    return gStats;
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
    // FIX: print the actual min/max/avg values, not overall_gc_content
    cout << "  > Minimum read: " << stats.minimum << " bp" << endl;
    cout << "  > Maximum read: " << stats.maximum << " bp" << endl;
    cout << "  > Avg length read: " << stats.avg_read_len << " bp" << endl;
    cout << "============================================" << endl;
}

int fastq_parser(const string& input_file, const string& output_file, const int qmin) {

    // Step 1: READING
    vector<FastqRecord> data = parseFastqFile(input_file, output_file, qmin);

    if (data.empty()) {
        cout << "No sequences found." << endl;
        return 0;
    }

    // Step 2: CALCULATING
    FastqGlobalStats stats = computeStatistics(data);

    // Step 3: PRINTING STATS
    printStatistics(stats);

    return 0;
}

//must check that,
// each record is complete,
// sequence and quality lines have equal length
// must compute minimum, maximum and average read length

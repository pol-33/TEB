#include "fastq_parser.hpp"

using namespace std;

// Parsing of the Fastq file
vector<FastqRecord> parseFastqFile(const string& infile, const string& outfile, const int qmin = 0) {
    vector<FastqRecord> records;
    ifstream in(infile);
    ofstream out(outfile);

    if (!in.is_open()) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return records;
    } else if (!out.is_open()) {
        cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
        return records;
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
                records.push_back(rec);
            }

            current_header = line.substr(1);
            current_sequence = "";
            out << current_header << "\n";
        } else {
            current_sequence += line;
            getline(in, line); //"+\n" line
            getline(in, line);
            sequence_quality += line;
            if (current_sequence.length() != sequence_quality.length()) {
                throw runtime_error("[parseFastqFile] sequence length != quality sequence length");
            }
            if (qmin > 0) {
                int last_good_pos = 0
                for (int i = current_sequence.length()-1; i >= 0; i--) {
                    int score = sequence_quality[i]-33;
                    if (score >= qmin)
                }
            }
        }
    }

    // We save the last record if it exists after the loop ends
    if (!current_header.empty()) {
        FastqRecord rec;
        rec.header = current_header;
        rec.sequence = current_sequence;
        records.push_back(rec);
    }

    in.close();
    out.close();
    return records;
}

// Compute statistics
GlobalStats computeStatistics(const vector<FastqRecord>& records) {
    GlobalStats gStats;
    gStats.num_sequences = records.size();
    gStats.total_length = 0;
    gStats.total_gc_count = 0;

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
    }

    // Calculate global GC
    gStats.overall_gc_content = (gStats.total_length > 0)
                                ? (double)gStats.total_gc_count / gStats.total_length * 100.0
                                : 0.0;

    return gStats;
}

// Print the statistics
void printStatistics(const GlobalStats& stats) {
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
    cout << "  > Minimum read: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Maximum read: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Avg length read: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "============================================" << endl;
}

int Fastq_parser(const string& input_file, const string& output_file) {

    // Step 1: READING
    vector<FastqRecord> data = parseFastqFile(input_file, output_file);

    if (data.empty()) {
        cout << "No sequences found." << endl;
        return 0;
    }

    // Step 2: CALCULATING
    GlobalStats stats = computeStatistics(data);

    // Step 3: PRINTING STATS
    printStatistics(stats);

    return 0;
}

//must check that, 
// each record is complete, 
// sequence and quality lines have equal length
// must compute minimum, maximum and average read length

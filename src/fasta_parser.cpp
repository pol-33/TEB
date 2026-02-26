#include "fasta_parser.hpp"

using namespace std;

// Parsing of the FASTA file
vector<FastaRecord> parseFastaFile(const string& infile, const string& outfile) {
    vector<FastaRecord> records;
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

    while (getline(in, line)) {
        if (line.empty()) continue;

        // Character '>', indicates new header
        if (line[0] == '>') {
            // save previous content
            if (!current_header.empty()) {
                FastaRecord rec;
                rec.header = current_header;
                rec.sequence = current_sequence;
                records.push_back(rec);
            }

            current_header = line.substr(1);
            current_sequence = "";
            out << current_header << "\n";
        } else {
            // Sequence line -> concat
            current_sequence += line;
            out << line << "\n";
        }
    }

    // We save the last record if it exists after the loop ends
    if (!current_header.empty()) {
        FastaRecord rec;
        rec.header = current_header;
        rec.sequence = current_sequence;
        records.push_back(rec);
    }

    in.close();
    out.close();
    return records;
}

// Compute statistics
GlobalStats computeStatistics(const vector<FastaRecord>& records) {
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
    cout << "============================================" << endl;
}

int fasta_parser(const string& input_file, const string& output_file) {

    // Step 1: READING
    vector<FastaRecord> data = parseFastaFile(input_file, output_file);

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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>

using namespace std;

// Estructura per guardar cada registre FASTA
struct FastaRecord {
    string header;
    string sequence;
};

struct SequenceStats {
    long long length;
    long long gc_count;
    double gc_content;
};

struct GlobalStats {
    int num_sequences;
    long long total_length;
    long long total_gc_count;
    double overall_gc_content;
    // Map: clau = ID de la seqüència, valor = les seves estadístiques
    map<string, SequenceStats> per_sequence; 
};


// Parsing of the FASTA file
vector<FastaRecord> parseFastaFile(const string& filename) {
    vector<FastaRecord> records;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: The file " << filename << " could not be opened." << endl;
        return records;
    }

    string line;
    string current_header = "";
    string current_sequence = "";

    while (getline(file, line)) {
        if (line.empty()) continue;

        // Si trobem el caràcter '>', se'ns està indicant que comença un nou registre
        if (line[0] == '>') {
            // Si ja teníem una capçalera llegida, vol dir que ha finalitzat la lectura del record anterior
            // Per tant, guardem el registre anterior
            if (!current_header.empty()) {
                FastaRecord rec;
                rec.header = current_header;
                rec.sequence = current_sequence;
                records.push_back(rec);
            }

            // Reiniciem per al nou registre
            current_header = line.substr(1);
            current_sequence = "";
        } else {
            // És una línia de seqüència, la concatenem
            current_sequence += line;
        }
    }

    // We save the last record if it exists after the loop ends
    if (!current_header.empty()) {
        FastaRecord rec;
        rec.header = current_header;
        rec.sequence = current_sequence;
        records.push_back(rec);
    }

    file.close();
    return records;
}


// Compute statistics
GlobalStats computeStatistics(const vector<FastaRecord>& records) {
    GlobalStats gStats;
    gStats.num_sequences = records.size();
    gStats.total_length = 0;
    gStats.total_gc_count = 0;

    // Calculem les estadístiques per a cada seqüència registrada
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

        // Calcular GC individual
        sStats.gc_content = (sStats.length > 0) 
                            ? (double)sStats.gc_count / sStats.length * 100.0 
                            : 0.0;

        // Guardar al map i actualitzar totals
        gStats.per_sequence[rec.header] = sStats;
        gStats.total_length += sStats.length;
        gStats.total_gc_count += sStats.gc_count;
    }

    // Calcular GC global
    gStats.overall_gc_content = (gStats.total_length > 0) 
                                ? (double)gStats.total_gc_count / gStats.total_length * 100.0 
                                : 0.0;

    return gStats;
}

// Print the statistics
void printStatistics(const GlobalStats& stats) {
    cout << "============================================" << endl;
    cout << "       INFORME D'ANÀLISI GENÒMICA" << endl;
    cout << "============================================" << endl;
    cout << "Nombre de seqüències processades: " << stats.num_sequences << endl;
    cout << "--------------------------------------------" << endl;

    // Iterem sobre el map per mostrar cada seqüència
    // it.first és la clau (header), it.second és l'objecte SequenceStats
    for (auto const& [header, sStats] : stats.per_sequence) {
        cout << "ID: " << header << endl;
        cout << "  > Longitud: " << sStats.length << " bp" << endl;
        cout << "  > Contingut GC: " << fixed << setprecision(2) << sStats.gc_content << "%" << endl;
        cout << endl;
    }

    cout << "--------------------------------------------" << endl;
    cout << "RESUM GLOBAL:" << endl;
    cout << "  > Longitud total: " << stats.total_length << " bp" << endl;
    cout << "  > Contingut GC total: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "============================================" << endl;
}

int main() {
    string path = "../chr1.fna";
    
    // Pas 1: Llegir
    vector<FastaRecord> data = parseFastaFile(path);
    
    if (data.empty()) {
        cout << "No s'han trobat seqüències." << endl;
        return 0;
    }

    // Pas 2: Calcular
    GlobalStats stats = computeStatistics(data);

    // Pas 3: Imprimir
    printStatistics(stats);

    return 0;
}

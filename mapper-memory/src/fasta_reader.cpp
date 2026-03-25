#include "fasta_reader.hpp"

#include <cctype>
#include <fstream>
#include <stdexcept>

#include "nucleotide.hpp"

namespace mapper_memory {

FastaData load_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA file: " + path);
    }

    FastaData data;
    std::string line;
    ChromInfo current;
    bool have_current = false;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            if (have_current) {
                data.chromosomes.push_back(current);
            }
            current = ChromInfo{};
            current.name = line.substr(1);
            current.offset = data.genome.size();
            current.length = 0;
            have_current = true;
            continue;
        }
        if (!have_current) {
            throw std::runtime_error("FASTA sequence found before any header in: " + path);
        }
        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
                continue;
            }
            data.genome.push_back(normalize_base(ch));
            ++current.length;
        }
    }

    if (have_current) {
        data.chromosomes.push_back(current);
    }
    if (data.chromosomes.empty()) {
        throw std::runtime_error("no chromosomes found in FASTA file: " + path);
    }

    return data;
}

}  // namespace mapper_memory

#include "fasta.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

std::string read_fasta_chrom_name(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("read_fasta_chrom_name: cannot open '" + path + "'");

    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (line[0] == '>') {
            std::string name = line.substr(1);
            const auto sp = name.find_first_of(" \t");
            if (sp != std::string::npos)
                name.resize(sp);
            return name;
        }
    }

    throw std::runtime_error("read_fasta_chrom_name: no FASTA header found in '" + path + "'");
}

// Read a (possibly multi-FASTA) file and return all sequence data concatenated
// into a single string.  Header lines (starting with '>') are skipped.
// Both LF and CRLF line endings are handled.
std::string read_fasta_sequence(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open())
        throw std::runtime_error("read_fasta_sequence: cannot open '" + path + "'");

    std::string result;
    std::string line;
    while (std::getline(in, line)) {
        // Strip Windows carriage return
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty() || line[0] == '>') continue;
        result += line;
    }

    std::cerr << "[fasta] Loaded " << result.size() << " bases\n";
    return result;
}

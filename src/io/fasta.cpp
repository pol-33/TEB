#include "fasta.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

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

#include "sam.hpp"

#include <fstream>
#include <stdexcept>
#include <string>

// SamWriter uses a private std::ofstream.  We add `out_` to sam.hpp (companion
// edit below).

SamWriter::SamWriter(const std::string& path) {
    out_.open(path);
    if (!out_.is_open())
        throw std::runtime_error("SamWriter: cannot open '" + path + "'");
}

void SamWriter::write_header(const std::string& ref_name, uint32_t ref_len) {
    out_ << "@HD\tVN:1.6\tSO:unsorted\n";
    out_ << "@SQ\tSN:" << ref_name
         << "\tLN:" << ref_len << "\n";
    if (!out_) throw std::runtime_error("SamWriter::write_header: write failed");
}

void SamWriter::write_alignment(const std::string& read_name,
                                const std::string& chrom,
                                uint32_t           pos,
                                const std::string& cigar,
                                const std::string& seq,
                                const std::string& qual,
                                const std::string& alt_chrom,
                                uint32_t           alt_pos,
                                const std::string& alt_cigar) {
    // Simplified alignment line format:
    // read_name chrom pos cigar seq qual [ALT:chrom,pos,cigar]
    out_ << read_name << '\t'
         << chrom << '\t'
         << pos << '\t'
         << cigar << '\t'
         << seq << '\t'
         << qual;

    if (!alt_chrom.empty() && alt_pos > 0 && !alt_cigar.empty()) {
        out_ << '\t'
             << "ALT:" << alt_chrom
             << ',' << alt_pos
             << ',' << alt_cigar;
    }
    out_ << '\n';

    if (!out_) throw std::runtime_error("SamWriter::write_alignment: write failed");
}

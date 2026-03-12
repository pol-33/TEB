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
                                uint32_t           pos,
                                const std::string& seq,
                                int                mismatches,
                                bool               mapped) {
    // FLAG
    const int flag = mapped ? 0 : 4;

    // Reference name and 1-based position
    const std::string rname = mapped ? "*" : "*";   // refined below
    const uint32_t    sam_pos = mapped ? pos + 1 : 0;
    const int         mapq    = mapped ? 60 : 0;
    const std::string cigar   = mapped
        ? std::to_string(seq.size()) + "M"
        : "*";
    const std::string nm_tag  = mapped
        ? "NM:i:" + std::to_string(mismatches)
        : "NM:i:0";

    // SAM columns: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS
    out_ << read_name      << '\t'   // QNAME
         << flag           << '\t'   // FLAG
         << "ref"          << '\t'   // RNAME  (single-ref index, always "ref")
         << sam_pos        << '\t'   // POS (1-based; 0 = unmapped)
         << mapq           << '\t'   // MAPQ
         << cigar          << '\t'   // CIGAR
         << "*"            << '\t'   // RNEXT
         << '0'            << '\t'   // PNEXT
         << '0'            << '\t'   // TLEN
         << seq            << '\t'   // SEQ
         << "*"            << '\t'   // QUAL
         << nm_tag         << '\n';  // optional NM tag
    if (!out_) throw std::runtime_error("SamWriter::write_alignment: write failed");
}

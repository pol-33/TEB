#pragma once

#include <cstdint>
#include <fstream>
#include <string>

class SamWriter {
public:
    explicit SamWriter(const std::string& path);
    void write_header(const std::string& ref_name, uint32_t ref_len);
    void write_alignment(
        const std::string& read_name,
        const std::string& chrom,
        uint32_t           pos,
        const std::string& cigar,
        const std::string& seq,
        const std::string& qual,
        const std::string& alt_chrom = "",
        uint32_t           alt_pos   = 0,
        const std::string& alt_cigar = ""
    );
private:
    std::ofstream out_;
};

#pragma once

#include <cstdint>
#include <string>

class SamWriter {
public:
    explicit SamWriter(const std::string& path);
    void write_header(const std::string& ref_name, uint32_t ref_len);
    void write_alignment(const std::string& read_name, uint32_t pos,
                         const std::string& seq, int mismatches, bool mapped);
};

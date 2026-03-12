#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

class GenomeStorage {
public:
    void encode(const std::string& seq);
    char get_base(uint32_t i) const;
    std::string substr(uint32_t pos, uint32_t len) const;
    bool is_n(uint32_t i) const;
    uint32_t size() const;
    const uint8_t*  packed_data() const;
    const uint64_t* n_mask_data() const;
    void write_binary(std::ostream& os) const;
    void read_binary(std::istream& is);

private:
    uint32_t              size_   = 0;
    std::vector<uint8_t>  packed_;   // 2-bit packed bases, 4 per byte, LSB-first
    std::vector<uint64_t> n_mask_;   // one bit per position; 1 = N
};

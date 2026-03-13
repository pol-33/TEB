#include "genome.hpp"
#include "binary.hpp"
#include "../util/dna_encoding.hpp"

#include <istream>
#include <ostream>
#include <stdexcept>

/*
 * Encoding: A=00  C=01  G=10  T=11  (2 bits each, 4 bases per byte, LSB-first)
 * N → stored as A (00) in packed_, but bit set in n_mask_.
 *
 * Inline test — encode("ACGTNACGT"), 9 bases:
 *
 *   pos:   0    1    2    3    4(N) 5    6    7    8
 *   2-bit: 00   01   10   11   00   00   01   10   11
 *
 *   Byte 0 (pos 0-3):  bits[7:6]=T=11  [5:4]=G=10  [3:2]=C=01  [1:0]=A=00
 *                    = 0b11_10_01_00 = 0xE4
 *   Byte 1 (pos 4-7):  bits[7:6]=G=10  [5:4]=C=01  [3:2]=A=00  [1:0]=N→A=00
 *                    = 0b10_01_00_00 = 0x90
 *   Byte 2 (pos 8):    bits[1:0]=T=11, upper bits unused (0)
 *                    = 0b00_00_00_11 = 0x03
 *
 *   packed_  = { 0xE4, 0x90, 0x03 }
 *   n_mask_[0] = 1 << 4 = 0x0000000000000010  (only position 4 is N)
 */

// byte value → 4 decoded chars (LSB-first sub-position order).
static constexpr auto make_lut() noexcept {
    struct LutTable { char data[256][4]; };
    LutTable t{};
    const char ALPHA[4] = {'A', 'C', 'G', 'T'};
    for (int v = 0; v < 256; ++v) {
        t.data[v][0] = ALPHA[(v >> 0) & 3];
        t.data[v][1] = ALPHA[(v >> 2) & 3];
        t.data[v][2] = ALPHA[(v >> 4) & 3];
        t.data[v][3] = ALPHA[(v >> 6) & 3];
    }
    return t;
}
static constexpr auto LUT = make_lut();

// ---------- public API ---------------------------------------------------- //

void GenomeStorage::encode(const std::string& seq) {
    const uint32_t len = static_cast<uint32_t>(seq.size());
    size_ = len;
    packed_.assign((len + 3) / 4, uint8_t{0});
    n_mask_.assign((len + 63) / 64, uint64_t{0});

    for (uint32_t i = 0; i < len; ++i) {
        const auto c = static_cast<unsigned char>(seq[i]);
        if (c == 'N' || c == 'n')
            n_mask_[i / 64] |= (uint64_t{1} << (i % 64));
        packed_[i / 4] |= static_cast<uint8_t>(dna::ENC_2BIT.v[c] << ((i % 4) * 2));
    }
}

char GenomeStorage::get_base(uint32_t i) const {
    if (i >= size_) throw std::out_of_range("GenomeStorage::get_base: index out of range");
    if ((n_mask_[i / 64] >> (i % 64)) & 1) return 'N';
    return LUT.data[packed_[i / 4]][i % 4];
}

std::string GenomeStorage::substr(uint32_t pos, uint32_t len) const {
    if (pos + len > size_) throw std::out_of_range("GenomeStorage::substr: range out of bounds");
    std::string out(len, '\0');
    for (uint32_t i = 0; i < len; ++i)
        out[i] = get_base(pos + i);
    return out;
}

bool GenomeStorage::is_n(uint32_t i) const {
    if (i >= size_) throw std::out_of_range("GenomeStorage::is_n: index out of range");
    return ((n_mask_[i / 64] >> (i % 64)) & 1) != 0;
}

uint32_t GenomeStorage::size() const { return size_; }

const uint8_t* GenomeStorage::packed_data() const {
    return packed_.empty() ? nullptr : packed_.data();
}

const uint64_t* GenomeStorage::n_mask_data() const {
    return n_mask_.empty() ? nullptr : n_mask_.data();
}

// Binary format: [size : uint64_t][packed : vec<uint8_t>][n_mask : vec<uint64_t>]
void GenomeStorage::write_binary(std::ostream& os) const {
    bio::write_pod(os, static_cast<uint64_t>(size_));
    bio::write_vec(os, packed_);
    bio::write_vec(os, n_mask_);
}

void GenomeStorage::read_binary(std::istream& is) {
    uint64_t len = 0;
    bio::read_pod(is, len);
    size_ = static_cast<uint32_t>(len);
    bio::read_vec(is, packed_);
    bio::read_vec(is, n_mask_);
}

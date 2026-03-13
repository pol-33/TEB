#pragma once

#include <cstdint>

namespace dna {

struct EncodeTable2Bit {
    uint8_t v[256]{};
};

constexpr EncodeTable2Bit make_encode_table_2bit() noexcept {
    EncodeTable2Bit t{};
    t.v[static_cast<unsigned>('A')] = t.v[static_cast<unsigned>('a')] = 0;
    t.v[static_cast<unsigned>('C')] = t.v[static_cast<unsigned>('c')] = 1;
    t.v[static_cast<unsigned>('G')] = t.v[static_cast<unsigned>('g')] = 2;
    t.v[static_cast<unsigned>('T')] = t.v[static_cast<unsigned>('t')] = 3;
    return t;
}

inline constexpr EncodeTable2Bit ENC_2BIT = make_encode_table_2bit();

} // namespace dna

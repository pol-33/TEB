#pragma once

#include <cstdint>

int count_mismatches(
    const uint8_t*  genome_packed,
    const uint64_t* genome_n_mask,
    uint32_t        genome_pos,
    const uint8_t*  read_packed,
    uint32_t        read_len
);

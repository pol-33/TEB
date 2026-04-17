#include "reference.hpp"

#include <fstream>
#include <stdexcept>

#include "common.hpp"
#include "nucleotide.hpp"

namespace mapper_speed {

ReferenceData load_reference(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA: " + path);
    }

    ReferenceData data;
    data.checksum = fnv1a_init();

    std::string line;
    ChromosomeRecord current;
    bool in_sequence = false;

    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            if (in_sequence) {
                data.chromosomes.push_back(current);
            }
            current = ChromosomeRecord{};
            const std::size_t cut = line.find_first_of(" \t", 1u);
            current.name = line.substr(1u, cut == std::string::npos ? std::string::npos : cut - 1u);
            current.start = data.genome_length;
            current.length = 0;
            in_sequence = true;
            continue;
        }
        if (!in_sequence) {
            throw std::runtime_error("FASTA sequence encountered before header");
        }
        for (char ch : line) {
            const char normalized = normalize_base(ch);
            const uint8_t code = base_to_code(normalized);
            packed_set(data.packed_bases, data.genome_length, code == kBaseCodeInvalid ? 0u : code);
            if (normalized == 'N') {
                bitset_set(data.n_mask_words, data.genome_length);
            }
            data.checksum = fnv1a_extend(data.checksum, &normalized, 1u);
            ++data.genome_length;
            ++current.length;
        }
    }

    if (in_sequence) {
        data.chromosomes.push_back(current);
    }
    if (data.chromosomes.empty()) {
        throw std::runtime_error("no chromosomes found in FASTA");
    }
    data.n_mask_words.resize((static_cast<std::size_t>(data.genome_length) + 63u) / 64u, 0);
    return data;
}

}  // namespace mapper_speed

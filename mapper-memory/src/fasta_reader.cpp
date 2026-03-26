#include "fasta_reader.hpp"

#include <cctype>
#include <stdexcept>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

std::string normalize_name(const std::string& raw_name) {
    const std::size_t cut = raw_name.find_first_of(" \t");
    return raw_name.substr(0, cut);
}

void append_normalized_sequence(std::vector<uint8_t>& packed_genome,
                                std::vector<uint8_t>& text,
                                uint64_t& genome_length,
                                uint64_t& chrom_length,
                                const std::string& line) {
    for (char ch : line) {
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            continue;
        }
        const uint8_t rank = rank_from_base(ch);
        append_packed_rank(packed_genome, genome_length, rank);
        text.push_back(rank);
        ++genome_length;
        ++chrom_length;
    }
}

}  // namespace

FastaReader::FastaReader(const std::string& path) : in_(path) {
    if (!in_) {
        throw std::runtime_error("failed to open FASTA file: " + path);
    }
}

bool FastaReader::next(FastaChromosome& chrom) {
    chrom = FastaChromosome{};
    std::string line;

    if (pending_header_.empty()) {
        while (std::getline(in_, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (!line.empty() && line.front() == '>') {
                pending_header_ = line.substr(1);
                break;
            }
        }
    }

    if (pending_header_.empty()) {
        return false;
    }

    chrom.name = normalize_name(pending_header_);
    chrom.offset = next_offset_;
    pending_header_.clear();

    while (std::getline(in_, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            pending_header_ = line.substr(1);
            break;
        }
        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) == 0) {
                chrom.sequence.push_back(normalize_base(ch));
            }
        }
    }

    if (chrom.sequence.empty()) {
        throw std::runtime_error("empty chromosome encountered in FASTA");
    }

    next_offset_ += chrom.sequence.size();
    return true;
}

FastaData load_fasta(const std::string& path) {
    FastaReader reader(path);
    FastaData data;
    FastaChromosome chrom;

    while (reader.next(chrom)) {
        ChromInfo info;
        info.name = chrom.name;
        info.offset = data.genome_length;
        info.length = chrom.sequence.size();
        data.chromosomes.push_back(info);

        uint64_t chrom_length = 0;
        append_normalized_sequence(data.packed_genome, data.text, data.genome_length, chrom_length, chrom.sequence);
        data.text.push_back(kSeparatorRank);
    }

    if (data.chromosomes.empty()) {
        throw std::runtime_error("no chromosomes found in FASTA file: " + path);
    }

    data.text.back() = kSentinelRank;
    return data;
}

}  // namespace mapper_memory

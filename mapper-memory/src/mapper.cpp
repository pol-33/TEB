#include <algorithm>
#include <array>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "alignment.hpp"
#include "backward_search.hpp"
#include "fasta_reader.hpp"
#include "fastq_reader.hpp"
#include "fm_index.hpp"
#include "nucleotide.hpp"

namespace {

struct Config {
    std::string index_path;
    std::string reference_path;
    std::string reads_path;
    std::string output_path;
    int max_errors = 0;
};

struct Candidate {
    uint64_t genome_pos = 0;
    int search_edit = 0;
};

struct VerifiedAlignment {
    mapper_memory::Alignment alignment;
    uint64_t genome_pos = 0;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-I" && i + 1 < argc) {
            config.index_path = argv[++i];
        } else if (arg == "-R" && i + 1 < argc) {
            config.reference_path = argv[++i];
        } else if (arg == "-i" && i + 1 < argc) {
            config.reads_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            config.output_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            config.max_errors = std::stoi(argv[++i]);
        } else {
            throw std::runtime_error("usage: mapper (-I genome.idx | -R genome.fa) -i reads.fastq -o output.sam -k <0..3>");
        }
    }

    if ((config.index_path.empty() && config.reference_path.empty()) ||
        (!config.index_path.empty() && !config.reference_path.empty()) ||
        config.reads_path.empty() || config.output_path.empty() ||
        config.max_errors < 0 || config.max_errors > 3) {
        throw std::runtime_error("usage: mapper (-I genome.idx | -R genome.fa) -i reads.fastq -o output.sam -k <0..3>");
    }

    return config;
}

bool better_alignment(const mapper_memory::Alignment& lhs, const mapper_memory::Alignment& rhs) {
    if (lhs.edit_dist != rhs.edit_dist) {
        return lhs.edit_dist < rhs.edit_dist;
    }
    if (lhs.chrom != rhs.chrom) {
        return lhs.chrom < rhs.chrom;
    }
    if (lhs.ref_pos != rhs.ref_pos) {
        return lhs.ref_pos < rhs.ref_pos;
    }
    return lhs.cigar < rhs.cigar;
}

void maybe_push_alignment(std::vector<VerifiedAlignment>& best, const VerifiedAlignment& candidate) {
    for (const VerifiedAlignment& existing : best) {
        if (existing.genome_pos == candidate.genome_pos) {
            return;
        }
    }
    best.push_back(candidate);
    std::sort(best.begin(), best.end(), [](const VerifiedAlignment& lhs, const VerifiedAlignment& rhs) {
        return better_alignment(lhs.alignment, rhs.alignment);
    });
    if (best.size() > 2U) {
        best.resize(2U);
    }
}

mapper_memory::Alignment verify_index_candidate(const mapper_memory::FMIndexView& index,
                                                const std::string& normalized_read,
                                                uint64_t genome_pos,
                                                int max_errors,
                                                mapper_memory::AlignmentWorkspace& workspace,
                                                std::string& ref_buffer) {
    mapper_memory::Alignment best;
    best.edit_dist = std::numeric_limits<int>::max();

    const uint64_t chrom_index = index.chromosome_index(genome_pos);
    const uint64_t chrom_end = index.chromosome_offset(chrom_index) + index.chromosome_length(chrom_index);
    const uint64_t min_ref_len = normalized_read.size() > static_cast<std::size_t>(max_errors)
        ? normalized_read.size() - static_cast<std::size_t>(max_errors)
        : 1U;
    const uint64_t max_ref_len = std::min<uint64_t>(chrom_end - genome_pos,
                                                    normalized_read.size() + static_cast<std::size_t>(max_errors));

    if (min_ref_len > max_ref_len) {
        return best;
    }

    ref_buffer.reserve(normalized_read.size() + static_cast<std::size_t>(max_errors));
    for (uint64_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        index.extract_reference(genome_pos, ref_len, ref_buffer);
        mapper_memory::Alignment current = mapper_memory::band_align(normalized_read, ref_buffer, max_errors, workspace);
        if (current.edit_dist <= max_errors && current.edit_dist < best.edit_dist) {
            current.chrom = std::string(index.chromosome_name(chrom_index));
            current.ref_pos = genome_pos + 1ULL;
            best = current;
        }
    }

    return best;
}

mapper_memory::Alignment verify_fallback_candidate(const mapper_memory::FastaData& fasta,
                                                   std::size_t chrom_index,
                                                   const std::string& normalized_read,
                                                   uint64_t genome_pos,
                                                   int max_errors,
                                                   mapper_memory::AlignmentWorkspace& workspace,
                                                   std::string& ref_buffer) {
    mapper_memory::Alignment best;
    best.edit_dist = std::numeric_limits<int>::max();

    const mapper_memory::ChromInfo& chrom = fasta.chromosomes[chrom_index];
    const uint64_t chrom_end = chrom.offset + chrom.length;
    const uint64_t min_ref_len = normalized_read.size() > static_cast<std::size_t>(max_errors)
        ? normalized_read.size() - static_cast<std::size_t>(max_errors)
        : 1U;
    const uint64_t max_ref_len = std::min<uint64_t>(chrom_end - genome_pos,
                                                    normalized_read.size() + static_cast<std::size_t>(max_errors));
    if (min_ref_len > max_ref_len) {
        return best;
    }

    ref_buffer.reserve(normalized_read.size() + static_cast<std::size_t>(max_errors));
    for (uint64_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        ref_buffer.assign(fasta.genome.data() + genome_pos, static_cast<std::size_t>(ref_len));
        mapper_memory::Alignment current = mapper_memory::band_align(normalized_read, ref_buffer, max_errors, workspace);
        if (current.edit_dist <= max_errors && current.edit_dist < best.edit_dist) {
            current.chrom = chrom.name;
            current.ref_pos = genome_pos + 1ULL;
            best = current;
        }
    }

    return best;
}

void collect_candidates(const mapper_memory::FMIndexView& index,
                        const std::vector<mapper_memory::SearchResult>& ranges,
                        std::size_t read_len,
                        int max_errors,
                        std::vector<Candidate>& candidates) {
    candidates.clear();
    std::unordered_set<uint64_t> seen;
    seen.reserve(64);

    for (const mapper_memory::SearchResult& result : ranges) {
        std::size_t emitted = 0;
        for (uint64_t row = result.sa_lo; row < result.sa_hi && candidates.size() < 32U && emitted < 16U; ++row) {
            const uint64_t sa_value = index.sa_value(row);
            if (sa_value >= index.genome_length()) {
                continue;
            }
            if (!index.same_chromosome_window(sa_value, read_len + static_cast<std::size_t>(max_errors))) {
                continue;
            }
            if (seen.insert(sa_value).second) {
                candidates.push_back(Candidate{sa_value, result.edit_dist});
                ++emitted;
            }
        }
    }

    std::sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
        if (lhs.search_edit != rhs.search_edit) {
            return lhs.search_edit < rhs.search_edit;
        }
        return lhs.genome_pos < rhs.genome_pos;
    });
}

void write_record(std::ofstream& out,
                  const mapper_memory::Read& read,
                  const std::vector<VerifiedAlignment>& alignments) {
    if (alignments.empty()) {
        out << read.name << " * 0 * " << read.seq << ' ' << read.qual << '\n';
        return;
    }

    const mapper_memory::Alignment& primary = alignments.front().alignment;
    out << read.name << ' '
        << primary.chrom << ' '
        << primary.ref_pos << ' '
        << primary.cigar << ' '
        << read.seq << ' '
        << read.qual;
    if (alignments.size() > 1U) {
        const mapper_memory::Alignment& alt = alignments[1].alignment;
        out << " ALT:" << alt.chrom << ',' << alt.ref_pos << ',' << alt.cigar;
    }
    out << '\n';
}

void run_index_mode(const Config& config) {
    mapper_memory::FMIndexView index(config.index_path);
    mapper_memory::FastqReader reader(config.reads_path);
    std::ofstream out(config.output_path);
    if (!out) {
        throw std::runtime_error("failed to open SAM output: " + config.output_path);
    }
    std::vector<char> output_buffer(4 * 1024 * 1024);
    out.rdbuf()->pubsetbuf(output_buffer.data(), static_cast<std::streamsize>(output_buffer.size()));

    mapper_memory::AlignmentWorkspace workspace;
    mapper_memory::Read read;
    std::vector<mapper_memory::SearchResult> ranges;
    ranges.reserve(32);
    std::vector<Candidate> candidates;
    candidates.reserve(32);
    std::vector<VerifiedAlignment> best_alignments;
    best_alignments.reserve(2);
    std::string ref_buffer;

    while (reader.next(read)) {
        const std::string normalized = mapper_memory::normalize_sequence(read.seq);
        mapper_memory::inexact_search(index, normalized, config.max_errors, ranges);
        collect_candidates(index, ranges, normalized.size(), config.max_errors, candidates);

        best_alignments.clear();
        for (const Candidate& candidate : candidates) {
            mapper_memory::Alignment alignment =
                verify_index_candidate(index, normalized, candidate.genome_pos, config.max_errors, workspace, ref_buffer);
            if (alignment.edit_dist <= config.max_errors) {
                maybe_push_alignment(best_alignments, VerifiedAlignment{alignment, candidate.genome_pos});
            }
        }

        write_record(out, read, best_alignments);
    }
}

void run_fallback_mode(const Config& config) {
    std::cerr << "[mapper] direct -R mode uses a correctness-first genome scan and is not competition-optimized\n";

    mapper_memory::FastaData fasta = mapper_memory::load_fasta(config.reference_path);
    mapper_memory::FastqReader reader(config.reads_path);
    std::ofstream out(config.output_path);
    if (!out) {
        throw std::runtime_error("failed to open SAM output: " + config.output_path);
    }
    std::vector<char> output_buffer(4 * 1024 * 1024);
    out.rdbuf()->pubsetbuf(output_buffer.data(), static_cast<std::streamsize>(output_buffer.size()));

    mapper_memory::AlignmentWorkspace workspace;
    mapper_memory::Read read;
    std::vector<VerifiedAlignment> best_alignments;
    best_alignments.reserve(2);
    std::string ref_buffer;

    while (reader.next(read)) {
        const std::string normalized = mapper_memory::normalize_sequence(read.seq);
        best_alignments.clear();

        for (std::size_t chrom_index = 0; chrom_index < fasta.chromosomes.size(); ++chrom_index) {
            const mapper_memory::ChromInfo& chrom = fasta.chromosomes[chrom_index];
            if (chrom.length == 0) {
                continue;
            }
            for (uint64_t offset = 0; offset < chrom.length; ++offset) {
                const uint64_t genome_pos = chrom.offset + offset;
                mapper_memory::Alignment alignment =
                    verify_fallback_candidate(fasta, chrom_index, normalized, genome_pos, config.max_errors, workspace, ref_buffer);
                if (alignment.edit_dist <= config.max_errors) {
                    maybe_push_alignment(best_alignments, VerifiedAlignment{alignment, genome_pos});
                }
            }
        }

        write_record(out, read, best_alignments);
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::ios::sync_with_stdio(false);
        const Config config = parse_args(argc, argv);
        if (!config.index_path.empty()) {
            run_index_mode(config);
        } else {
            run_fallback_mode(config);
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mapper error: " << e.what() << '\n';
        return 1;
    }
}

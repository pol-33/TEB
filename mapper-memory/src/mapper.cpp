#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "alignment.hpp"
#include "fastq_reader.hpp"
#include "fm_index.hpp"
#include "fm_search.hpp"
#include "memory_stats.hpp"
#include "nucleotide.hpp"

namespace {

struct Config {
    std::string index_path;
    std::string reads_path;
    std::string output_path;
    int max_errors = 0;
};

struct VerifiedAlignment {
    mapper_memory::Alignment alignment;
    std::size_t chrom_index = 0;
    uint64_t text_pos = 0;
};

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-I" && i + 1 < argc) {
            config.index_path = argv[++i];
        } else if (arg == "-i" && i + 1 < argc) {
            config.reads_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            config.output_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            config.max_errors = std::stoi(argv[++i]);
        } else if (arg == "-R") {
            throw std::runtime_error("direct -R mode is disabled for the contest-scale mapper; build an FM-index first");
        } else {
            throw std::runtime_error("usage: mapper -I genome.idx -i reads.fastq -o output.sam -k <0..3>");
        }
    }
    if (config.index_path.empty() || config.reads_path.empty() || config.output_path.empty() ||
        config.max_errors < 0 || config.max_errors > 3) {
        throw std::runtime_error("usage: mapper -I genome.idx -i reads.fastq -o output.sam -k <0..3>");
    }
    return config;
}

bool better_alignment(const VerifiedAlignment& lhs, const VerifiedAlignment& rhs) {
    if (lhs.alignment.edit_dist != rhs.alignment.edit_dist) {
        return lhs.alignment.edit_dist < rhs.alignment.edit_dist;
    }
    if (lhs.alignment.chrom != rhs.alignment.chrom) {
        return lhs.alignment.chrom < rhs.alignment.chrom;
    }
    if (lhs.alignment.ref_pos != rhs.alignment.ref_pos) {
        return lhs.alignment.ref_pos < rhs.alignment.ref_pos;
    }
    return lhs.alignment.cigar < rhs.alignment.cigar;
}

bool same_alignment_bucket(const VerifiedAlignment& lhs,
                           const VerifiedAlignment& rhs,
                           std::size_t read_length) {
    if (lhs.alignment.chrom != rhs.alignment.chrom) {
        return false;
    }
    const uint64_t lhs_pos = lhs.alignment.ref_pos;
    const uint64_t rhs_pos = rhs.alignment.ref_pos;
    const uint64_t delta = (lhs_pos > rhs_pos) ? (lhs_pos - rhs_pos) : (rhs_pos - lhs_pos);
    return delta <= static_cast<uint64_t>(read_length);
}

void maybe_push_alignment(std::vector<VerifiedAlignment>& best,
                          const VerifiedAlignment& candidate,
                          std::size_t read_length) {
    for (VerifiedAlignment& existing : best) {
        if ((existing.chrom_index == candidate.chrom_index && existing.text_pos == candidate.text_pos) ||
            same_alignment_bucket(existing, candidate, read_length)) {
            if (better_alignment(candidate, existing)) {
                existing = candidate;
                std::sort(best.begin(), best.end(), better_alignment);
            }
            return;
        }
    }
    best.push_back(candidate);
    std::sort(best.begin(), best.end(), better_alignment);
    if (best.size() > 2U) {
        best.resize(2U);
    }
}

char complement_base(char base) {
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default: return 'A';
    }
}

std::string reverse_complement_sequence(const std::string& normalized) {
    std::string out(normalized.size(), 'A');
    for (std::size_t i = 0; i < normalized.size(); ++i) {
        out[i] = complement_base(normalized[normalized.size() - 1U - i]);
    }
    return out;
}

mapper_memory::Alignment verify_candidate(const mapper_memory::FMIndexView::ChromosomeView& chrom,
                                          const std::string& normalized_read,
                                          uint64_t text_pos,
                                          int max_errors,
                                          mapper_memory::AlignmentWorkspace& workspace,
                                          std::string& ref_buffer) {
    mapper_memory::Alignment best;
    best.edit_dist = std::numeric_limits<int>::max();

    const uint64_t min_ref_len = normalized_read.size() > static_cast<std::size_t>(max_errors)
        ? normalized_read.size() - static_cast<std::size_t>(max_errors)
        : 1U;
    const uint64_t max_ref_len = std::min<uint64_t>(chrom.length - text_pos,
                                                    normalized_read.size() + static_cast<std::size_t>(max_errors));
    if (min_ref_len > max_ref_len) {
        return best;
    }

    ref_buffer.reserve(normalized_read.size() + static_cast<std::size_t>(max_errors));
    for (uint64_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        chrom.extract_reference(text_pos, ref_len, ref_buffer);
        mapper_memory::Alignment current = mapper_memory::band_align(normalized_read, ref_buffer, max_errors, workspace);
        if (current.edit_dist <= max_errors && current.edit_dist < best.edit_dist) {
            current.chrom = std::string(chrom.name);
            current.ref_pos = text_pos + 1ULL;
            best = current;
        }
    }
    return best;
}

void evaluate_interval(const mapper_memory::FMIndexView::ChromosomeView& chrom,
                       std::size_t chrom_index,
                       const mapper_memory::SearchResult& result,
                       const std::string& normalized_read,
                       int max_errors,
                       mapper_memory::AlignmentWorkspace& workspace,
                       std::string& ref_buffer,
                       std::unordered_set<uint64_t>& seen_positions,
                       std::vector<VerifiedAlignment>& best_alignments) {
    const uint64_t max_row = std::min<uint64_t>(result.sa_hi, result.sa_lo + 4ULL);
    for (uint64_t row = result.sa_lo; row < max_row; ++row) {
        const uint64_t text_pos = chrom.locate(row);
        if (text_pos >= chrom.length) {
            continue;
        }
        if (!seen_positions.insert((static_cast<uint64_t>(chrom_index) << 32U) ^ text_pos).second) {
            continue;
        }
        if (!chrom.same_chromosome_window(text_pos, normalized_read.size())) {
            continue;
        }
        mapper_memory::Alignment alignment =
            verify_candidate(chrom, normalized_read, text_pos, max_errors, workspace, ref_buffer);
        if (alignment.edit_dist <= max_errors) {
            maybe_push_alignment(best_alignments, VerifiedAlignment{alignment, chrom_index, text_pos}, normalized_read.size());
        }
    }
}

void search_orientation(const mapper_memory::FMIndexView& index,
                        const std::string& normalized_read,
                        int max_errors,
                        mapper_memory::AlignmentWorkspace& workspace,
                        std::string& ref_buffer,
                        std::unordered_set<uint64_t>& seen_positions,
                        std::vector<mapper_memory::SearchResult>& search_results,
                        std::vector<VerifiedAlignment>& best_alignments) {
    for (std::size_t chrom_index = 0; chrom_index < index.chromosome_count(); ++chrom_index) {
        const auto& chrom = index.chromosome(chrom_index);
        search_results.clear();

        const mapper_memory::SAInterval exact = mapper_memory::exact_search(chrom, normalized_read);
        if (exact.lo < exact.hi) {
            search_results.push_back(mapper_memory::SearchResult{exact.lo, exact.hi, 0});
        } else if (max_errors > 0) {
            mapper_memory::inexact_search(chrom, normalized_read, max_errors, search_results, 32U);
        }
        if (search_results.empty()) {
            continue;
        }

        const int best_interval_edit = search_results.front().edit_dist;
        for (const mapper_memory::SearchResult& result : search_results) {
            if (result.edit_dist > best_interval_edit) {
                break;
            }
            evaluate_interval(chrom, chrom_index, result, normalized_read, max_errors, workspace,
                              ref_buffer, seen_positions, best_alignments);
            if (best_alignments.size() == 2U && best_alignments.back().alignment.edit_dist == 0) {
                break;
            }
        }
    }
}

void write_record(std::ofstream& out,
                  const mapper_memory::Read& read,
                  const std::vector<VerifiedAlignment>& alignments) {
    if (alignments.empty()) {
        out << read.name << " * 0 * " << read.seq << ' ' << read.qual << '\n';
        return;
    }
    const mapper_memory::Alignment& primary = alignments.front().alignment;
    out << read.name << ' ' << primary.chrom << ' ' << primary.ref_pos << ' ' << primary.cigar
        << ' ' << read.seq << ' ' << read.qual;
    if (alignments.size() > 1U) {
        const mapper_memory::Alignment& alt = alignments[1].alignment;
        out << " ALT:" << alt.chrom << ',' << alt.ref_pos << ',' << alt.cigar;
    }
    out << '\n';
}

void report_peak_rss() {
    const mapper_memory::MemoryStats stats = mapper_memory::read_memory_stats();
    std::cerr << "[mapper] current RSS: " << mapper_memory::bytes_to_mebibytes(stats.current_rss_bytes)
              << " MiB, peak RSS: " << mapper_memory::bytes_to_mebibytes(stats.peak_rss_bytes) << " MiB\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::ios::sync_with_stdio(false);
        const Config config = parse_args(argc, argv);

        mapper_memory::FMIndexView index(config.index_path);
        mapper_memory::FastqReader reader(config.reads_path);
        std::ofstream out(config.output_path);
        if (!out) {
            throw std::runtime_error("failed to open output file: " + config.output_path);
        }
        std::vector<char> output_buffer(4 * 1024 * 1024);
        out.rdbuf()->pubsetbuf(output_buffer.data(), static_cast<std::streamsize>(output_buffer.size()));

        mapper_memory::AlignmentWorkspace workspace;
        mapper_memory::Read read;
        std::vector<mapper_memory::SearchResult> search_results;
        search_results.reserve(64);
        std::vector<VerifiedAlignment> best_alignments;
        best_alignments.reserve(2);
        std::unordered_set<uint64_t> seen_positions;
        seen_positions.reserve(1024);
        std::string ref_buffer;
        std::string normalized;
        std::string reverse_complement;

        uint64_t processed = 0;
        while (reader.next(read)) {
            normalized = mapper_memory::normalize_sequence(read.seq);
            best_alignments.clear();
            seen_positions.clear();

            search_orientation(index, normalized, config.max_errors, workspace, ref_buffer,
                               seen_positions, search_results, best_alignments);

            reverse_complement = reverse_complement_sequence(normalized);
            if (reverse_complement != normalized || best_alignments.empty()) {
                search_orientation(index, reverse_complement, config.max_errors, workspace, ref_buffer,
                                   seen_positions, search_results, best_alignments);
            }

            write_record(out, read, best_alignments);
            ++processed;
            if (processed % 100000ULL == 0ULL) {
                std::cerr << "[mapper] processed " << processed << " reads\n";
            }
        }

        std::cerr << "[mapper] processed " << processed << " reads total\n";
        report_peak_rss();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mapper error: " << e.what() << '\n';
        return 1;
    }
}

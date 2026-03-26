#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "alignment.hpp"
#include "fastq_reader.hpp"
#include "nucleotide.hpp"
#include "seed_index.hpp"

namespace {

struct Config {
    std::string index_path;
    std::string reads_path;
    std::string output_path;
    int max_errors = 0;
};

struct SeedProbe {
    uint32_t code = 0;
    uint32_t read_offset = 0;
    uint32_t hit_count = 0;
};

struct Candidate {
    uint64_t genome_pos = 0;
    uint32_t supporting_hits = 0;
};

struct VerifiedAlignment {
    mapper_memory::Alignment alignment;
    uint64_t genome_pos = 0;
    uint32_t supporting_hits = 0;
};

constexpr std::size_t kMaxCandidates = 8192;
constexpr uint32_t kMaxSeedHitsToExpand = 250000;

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
            throw std::runtime_error("direct -R mode is disabled for the contest-scale mapper; build an index first with ./indexer -R genome.fa -I genome.idx");
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
    if (lhs.supporting_hits != rhs.supporting_hits) {
        return lhs.supporting_hits > rhs.supporting_hits;
    }
    if (lhs.alignment.chrom != rhs.alignment.chrom) {
        return lhs.alignment.chrom < rhs.alignment.chrom;
    }
    if (lhs.alignment.ref_pos != rhs.alignment.ref_pos) {
        return lhs.alignment.ref_pos < rhs.alignment.ref_pos;
    }
    return lhs.alignment.cigar < rhs.alignment.cigar;
}

void maybe_push_alignment(std::vector<VerifiedAlignment>& best, const VerifiedAlignment& candidate) {
    for (VerifiedAlignment& existing : best) {
        if (existing.genome_pos == candidate.genome_pos) {
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
        case 'A':
            return 'T';
        case 'C':
            return 'G';
        case 'G':
            return 'C';
        case 'T':
            return 'A';
        default:
            return 'A';
    }
}

std::string reverse_complement_sequence(const std::string& normalized) {
    std::string out;
    out.resize(normalized.size());
    for (std::size_t i = 0; i < normalized.size(); ++i) {
        out[i] = complement_base(normalized[normalized.size() - 1U - i]);
    }
    return out;
}

std::vector<std::pair<uint32_t, uint32_t>> make_segments(std::size_t read_length, int max_errors) {
    const std::size_t blocks = static_cast<std::size_t>(max_errors) + 1U;
    std::vector<std::pair<uint32_t, uint32_t>> segments;
    segments.reserve(blocks);

    const std::size_t base = read_length / blocks;
    const std::size_t extra = read_length % blocks;
    std::size_t cursor = 0;
    for (std::size_t i = 0; i < blocks; ++i) {
        const std::size_t len = base + (i < extra ? 1U : 0U);
        segments.emplace_back(static_cast<uint32_t>(cursor), static_cast<uint32_t>(len));
        cursor += len;
    }
    return segments;
}

bool choose_probe_for_segment(const mapper_memory::SeedIndexView& index,
                              const std::string& read,
                              uint32_t segment_offset,
                              uint32_t segment_length,
                              SeedProbe& best_probe) {
    if (segment_length < mapper_memory::kSeedLength) {
        return false;
    }

    bool found = false;
    best_probe.hit_count = std::numeric_limits<uint32_t>::max();
    for (uint32_t delta = 0; delta + mapper_memory::kSeedLength <= segment_length; ++delta) {
        const uint32_t read_offset = segment_offset + delta;
        const uint32_t code = mapper_memory::encode_seed(
            std::string_view(read.data() + read_offset, mapper_memory::kSeedLength));
        const uint32_t begin = index.bucket_begin(code);
        const uint32_t end = index.bucket_end(code);
        const uint32_t hits = end - begin;
        if (hits == 0) {
            continue;
        }
        if (!found || hits < best_probe.hit_count) {
            found = true;
            best_probe.code = code;
            best_probe.read_offset = read_offset;
            best_probe.hit_count = hits;
            if (hits == 1) {
                break;
            }
        }
    }

    return found;
}

void collect_probes(const mapper_memory::SeedIndexView& index,
                    const std::string& read,
                    int max_errors,
                    std::vector<SeedProbe>& probes) {
    probes.clear();
    const std::vector<std::pair<uint32_t, uint32_t>> segments = make_segments(read.size(), max_errors);
    for (const auto& segment : segments) {
        SeedProbe probe;
        if (choose_probe_for_segment(index, read, segment.first, segment.second, probe)) {
            probes.push_back(probe);
        }
    }
    std::sort(probes.begin(), probes.end(), [](const SeedProbe& lhs, const SeedProbe& rhs) {
        if (lhs.hit_count != rhs.hit_count) {
            return lhs.hit_count < rhs.hit_count;
        }
        return lhs.read_offset < rhs.read_offset;
    });
}

void collect_candidates(const mapper_memory::SeedIndexView& index,
                        const std::vector<SeedProbe>& probes,
                        std::size_t read_length,
                        int max_errors,
                        std::vector<Candidate>& candidates,
                        std::unordered_map<uint64_t, uint32_t>& support_counts) {
    candidates.clear();
    support_counts.clear();

    bool expanded_any = false;
    const uint64_t min_required_span =
        (read_length > static_cast<std::size_t>(max_errors)) ? read_length - static_cast<std::size_t>(max_errors) : 1U;
    for (const SeedProbe& probe : probes) {
        if (probe.hit_count > kMaxSeedHitsToExpand && expanded_any) {
            continue;
        }
        const uint32_t begin = index.bucket_begin(probe.code);
        const uint32_t end = index.bucket_end(probe.code);
        for (uint32_t i = begin; i < end && candidates.size() < kMaxCandidates; ++i) {
            const uint32_t seed_pos = index.position_at(i);
            for (int shift = -max_errors; shift <= max_errors && candidates.size() < kMaxCandidates; ++shift) {
                const int64_t start_signed =
                    static_cast<int64_t>(seed_pos) - static_cast<int64_t>(probe.read_offset) + static_cast<int64_t>(shift);
                if (start_signed < 0) {
                    continue;
                }
                const uint64_t start = static_cast<uint64_t>(start_signed);
                if (start + min_required_span > index.genome_length()) {
                    continue;
                }
                if (!index.same_chromosome_window(start, min_required_span)) {
                    continue;
                }
                auto [it, inserted] = support_counts.emplace(start, 1U);
                if (!inserted) {
                    ++it->second;
                    continue;
                }
                if (candidates.size() < kMaxCandidates) {
                    candidates.push_back(Candidate{start, 1U});
                }
            }
        }
        expanded_any = true;
        if (candidates.size() >= kMaxCandidates) {
            break;
        }
    }

    for (Candidate& candidate : candidates) {
        candidate.supporting_hits = support_counts[candidate.genome_pos];
    }
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
        if (lhs.supporting_hits != rhs.supporting_hits) {
            return lhs.supporting_hits > rhs.supporting_hits;
        }
        return lhs.genome_pos < rhs.genome_pos;
    });
}

void collect_candidates_by_scan(const mapper_memory::SeedIndexView& index,
                                std::size_t read_length,
                                int max_errors,
                                std::vector<Candidate>& candidates) {
    candidates.clear();
    const uint64_t min_required_span =
        (read_length > static_cast<std::size_t>(max_errors)) ? read_length - static_cast<std::size_t>(max_errors) : 1U;

    for (const auto& chrom : index.chromosomes()) {
        if (chrom.length < min_required_span) {
            continue;
        }
        const uint64_t last_start = chrom.offset + chrom.length - min_required_span;
        for (uint64_t pos = chrom.offset; pos <= last_start && candidates.size() < kMaxCandidates; ++pos) {
            candidates.push_back(Candidate{pos, 1U});
        }
        if (candidates.size() >= kMaxCandidates) {
            break;
        }
    }
}

mapper_memory::Alignment verify_candidate(const mapper_memory::SeedIndexView& index,
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
            current.ref_pos = genome_pos - index.chromosome_offset(chrom_index) + 1ULL;
            best = current;
        }
    }
    return best;
}

void search_orientation(const mapper_memory::SeedIndexView& index,
                        const std::string& normalized_read,
                        int max_errors,
                        std::vector<SeedProbe>& probes,
                        std::vector<Candidate>& candidates,
                        std::unordered_map<uint64_t, uint32_t>& support_counts,
                        mapper_memory::AlignmentWorkspace& workspace,
                        std::string& ref_buffer,
                        std::vector<VerifiedAlignment>& best_alignments) {
    collect_probes(index, normalized_read, max_errors, probes);
    if (!probes.empty()) {
        collect_candidates(index, probes, normalized_read.size(), max_errors, candidates, support_counts);
    } else {
        collect_candidates_by_scan(index, normalized_read.size(), max_errors, candidates);
    }

    for (const Candidate& candidate : candidates) {
        mapper_memory::Alignment alignment =
            verify_candidate(index, normalized_read, candidate.genome_pos, max_errors, workspace, ref_buffer);
        if (alignment.edit_dist <= max_errors) {
            maybe_push_alignment(best_alignments,
                                 VerifiedAlignment{alignment, candidate.genome_pos, candidate.supporting_hits});
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

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::ios::sync_with_stdio(false);
        const Config config = parse_args(argc, argv);

        mapper_memory::SeedIndexView index(config.index_path);
        mapper_memory::FastqReader reader(config.reads_path);
        std::ofstream out(config.output_path);
        if (!out) {
            throw std::runtime_error("failed to open output file: " + config.output_path);
        }
        std::vector<char> output_buffer(4 * 1024 * 1024);
        out.rdbuf()->pubsetbuf(output_buffer.data(), static_cast<std::streamsize>(output_buffer.size()));

        mapper_memory::AlignmentWorkspace workspace;
        mapper_memory::Read read;
        std::vector<SeedProbe> probes;
        probes.reserve(8);
        std::vector<Candidate> candidates;
        candidates.reserve(kMaxCandidates);
        std::unordered_map<uint64_t, uint32_t> support_counts;
        support_counts.reserve(kMaxCandidates * 2U);
        std::vector<VerifiedAlignment> best_alignments;
        best_alignments.reserve(2);
        std::string ref_buffer;
        std::string normalized;
        std::string reverse_complement;

        uint64_t processed = 0;
        while (reader.next(read)) {
            normalized = mapper_memory::normalize_sequence(read.seq);
            best_alignments.clear();
            search_orientation(index,
                               normalized,
                               config.max_errors,
                               probes,
                               candidates,
                               support_counts,
                               workspace,
                               ref_buffer,
                               best_alignments);

            reverse_complement = reverse_complement_sequence(normalized);
            if (reverse_complement != normalized || best_alignments.empty()) {
                search_orientation(index,
                                   reverse_complement,
                                   config.max_errors,
                                   probes,
                                   candidates,
                                   support_counts,
                                   workspace,
                                   ref_buffer,
                                   best_alignments);
            }

            write_record(out, read, best_alignments);
            ++processed;
            if (processed % 100000ULL == 0ULL) {
                std::cerr << "[mapper] processed " << processed << " reads\n";
            }
        }

        std::cerr << "[mapper] processed " << processed << " reads total\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mapper error: " << e.what() << '\n';
        return 1;
    }
}

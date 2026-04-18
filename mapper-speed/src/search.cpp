#include "search.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <string_view>

#include "common.hpp"
#include "nucleotide.hpp"

namespace mapper_speed {

namespace {

constexpr std::size_t kMaxDpFinalistsPerOrientation = 12;

uint32_t encode_seed_at(const std::string& read, uint32_t start, bool& valid) {
    uint32_t key = 0;
    valid = true;
    for (uint32_t i = 0; i < kSeedLength; ++i) {
        const uint8_t code = base_to_code(read[start + i]);
        if (code == kBaseCodeInvalid) {
            valid = false;
            return 0;
        }
        key = (key << 2u) | code;
    }
    return key;
}

std::size_t target_seed_count(int max_errors) {
    switch (max_errors) {
        case 0: return 4;
        case 1: return 6;
        case 2: return 7;
        default: return 8;
    }
}

}  // namespace

MapperEngine::MapperEngine(const IndexView& index)
    : index_(index),
      dispatch_(resolve_myers_dispatch()),
      candidate_table_(kCandidateTableSize) {
}

void MapperEngine::collect_seed_positions(std::size_t read_length,
                                          int max_errors,
                                          std::vector<uint32_t>& positions) const {
    positions.clear();
    if (read_length <= kSeedLength) {
        positions.push_back(0);
        return;
    }

    const uint32_t max_start = static_cast<uint32_t>(read_length - kSeedLength);
    const std::size_t target = std::min<std::size_t>(target_seed_count(max_errors), static_cast<std::size_t>(max_start) + 1u);
    const uint32_t partitions = static_cast<uint32_t>(std::max(1, max_errors + 1));

    positions.push_back(0);
    positions.push_back(max_start);
    for (uint32_t i = 0; i < partitions; ++i) {
        const uint32_t left = static_cast<uint32_t>((static_cast<uint64_t>(i) * max_start) / partitions);
        const uint32_t right = static_cast<uint32_t>((static_cast<uint64_t>(i + 1u) * max_start) / partitions);
        positions.push_back(left);
        positions.push_back((left + right) / 2u);
        positions.push_back(right);
    }

    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    if (positions.size() > target) {
        std::vector<uint32_t> trimmed;
        trimmed.reserve(target);
        for (std::size_t i = 0; i < target; ++i) {
            const std::size_t idx = (i * (positions.size() - 1u)) / (target - 1u);
            trimmed.push_back(positions[idx]);
        }
        std::sort(trimmed.begin(), trimmed.end());
        trimmed.erase(std::unique(trimmed.begin(), trimmed.end()), trimmed.end());
        positions.swap(trimmed);
    }

    if (positions.size() < target) {
        for (uint32_t pos = 0; pos <= max_start && positions.size() < target; ++pos) {
            if (!std::binary_search(positions.begin(), positions.end(), pos)) {
                positions.push_back(pos);
            }
        }
        std::sort(positions.begin(), positions.end());
    }
}

void MapperEngine::collect_seeds(const std::string& normalized_read,
                                 int max_errors,
                                 std::vector<SeedSpec>& seeds) {
    seeds.clear();
    scratch_seed_positions_.clear();
    collect_seed_positions(normalized_read.size(), max_errors, scratch_seed_positions_);

    bool any_under_cutoff = false;
    SeedSpec best_over_cutoff{};
    best_over_cutoff.frequency = std::numeric_limits<uint32_t>::max();

    for (uint32_t pos : scratch_seed_positions_) {
        bool valid = false;
        const uint32_t key = encode_seed_at(normalized_read, pos, valid);
        if (!valid) {
            continue;
        }
        const uint32_t freq = index_.occurrence_count(key);
        if (freq == 0) {
            continue;
        }
        SeedSpec spec{pos, key, freq};
        if (freq <= kHighFreqSkip) {
            any_under_cutoff = true;
            seeds.push_back(spec);
        } else if (freq <= kHighFreqAllowFallback && freq < best_over_cutoff.frequency) {
            best_over_cutoff = spec;
        }
    }

    if (!any_under_cutoff && best_over_cutoff.frequency != std::numeric_limits<uint32_t>::max()) {
        seeds.push_back(best_over_cutoff);
    }

    std::sort(seeds.begin(), seeds.end(), [](const SeedSpec& lhs, const SeedSpec& rhs) {
        if (lhs.frequency != rhs.frequency) {
            return lhs.frequency < rhs.frequency;
        }
        return lhs.read_offset < rhs.read_offset;
    });
}

void MapperEngine::generate_candidates(const std::vector<SeedSpec>& seeds,
                                       std::size_t read_length,
                                       int max_errors,
                                       std::vector<CandidateInfo>& starts) {
    starts.clear();
    for (CandidateSlot& slot : candidate_table_) {
        slot.occupied = false;
        slot.support = 0;
        slot.best_seed_freq = 0xFFFFFFFFu;
        slot.start = 0;
    }

    std::vector<uint32_t> per_seed;
    per_seed.reserve(4096);

    for (const SeedSpec& seed : seeds) {
        per_seed.clear();
        const auto range = index_.positions_for(seed.key);
        for (const uint32_t* it = range.first; it != range.second; ++it) {
            const int64_t ref_seed_pos = static_cast<int64_t>(*it);
            for (int delta = -max_errors; delta <= max_errors; ++delta) {
                const int64_t candidate = ref_seed_pos - static_cast<int64_t>(seed.read_offset) + delta;
                if (candidate < 0) {
                    continue;
                }
                const uint64_t start = static_cast<uint64_t>(candidate);
                const uint64_t min_ref_len = read_length > static_cast<std::size_t>(max_errors)
                    ? read_length - static_cast<std::size_t>(max_errors)
                    : 1u;
                if (start + min_ref_len > index_.genome_length()) {
                    continue;
                }
                per_seed.push_back(static_cast<uint32_t>(start));
            }
        }

        std::sort(per_seed.begin(), per_seed.end());
        per_seed.erase(std::unique(per_seed.begin(), per_seed.end()), per_seed.end());

        for (uint32_t start : per_seed) {
            std::size_t slot_idx = start & (candidate_table_.size() - 1u);
            while (true) {
                CandidateSlot& slot = candidate_table_[slot_idx];
                if (!slot.occupied) {
                    slot.occupied = true;
                    slot.start = start;
                    slot.support = 1;
                    slot.best_seed_freq = seed.frequency;
                    break;
                }
                if (slot.start == start) {
                    ++slot.support;
                    slot.best_seed_freq = std::min(slot.best_seed_freq, seed.frequency);
                    break;
                }
                slot_idx = (slot_idx + 1u) & (candidate_table_.size() - 1u);
            }
        }
    }

    std::vector<CandidateSlot> ranked;
    ranked.reserve(candidate_table_.size());
    for (const CandidateSlot& slot : candidate_table_) {
        if (slot.occupied) {
            ranked.push_back(slot);
        }
    }

    std::sort(ranked.begin(), ranked.end(), [](const CandidateSlot& lhs, const CandidateSlot& rhs) {
        if (lhs.support != rhs.support) {
            return lhs.support > rhs.support;
        }
        if (lhs.best_seed_freq != rhs.best_seed_freq) {
            return lhs.best_seed_freq < rhs.best_seed_freq;
        }
        return lhs.start < rhs.start;
    });

    const std::size_t limit = std::min<std::size_t>(kMaxVerifyPerOrientation, ranked.size());
    for (std::size_t i = 0; i < limit; ++i) {
        starts.push_back(CandidateInfo{
            ranked[i].start,
            ranked[i].best_seed_freq,
            ranked[i].support
        });
    }
}

int MapperEngine::prefilter_candidate(const std::string& oriented_read,
                                      const MyersQuery& query,
                                      const CandidateInfo& candidate,
                                      int max_errors) {
    const uint32_t global_start = candidate.start;
    const std::size_t chrom_index = index_.chromosome_for_position(global_start);
    const auto& chrom = index_.chromosome(chrom_index);
    const uint64_t local_pos = static_cast<uint64_t>(global_start) - chrom.start;
    const uint32_t min_ref_len = oriented_read.size() > static_cast<std::size_t>(max_errors)
        ? static_cast<uint32_t>(oriented_read.size() - static_cast<std::size_t>(max_errors))
        : 1u;
    const uint32_t max_ref_len = std::min<uint32_t>(
        static_cast<uint32_t>(oriented_read.size() + static_cast<std::size_t>(max_errors)),
        chrom.length - static_cast<uint32_t>(local_pos));

    if (min_ref_len > max_ref_len) {
        return std::numeric_limits<int>::max();
    }

    int best_score = std::numeric_limits<int>::max();
    for (uint32_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        index_.extract_sequence(global_start, ref_len, ref_buffer_);
        if (ref_len == oriented_read.size()) {
            const uint32_t mismatches = count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), ref_len);
            if (mismatches == 0u) {
                return 0;
            }
            best_score = std::min(best_score, static_cast<int>(mismatches));
            if (static_cast<int>(mismatches) <= max_errors) {
                continue;
            }
            const int score = bounded_edit_distance(dispatch_, query, ref_buffer_, max_errors);
            best_score = std::min(best_score, score);
        } else {
            const int score = banded_score_only(oriented_read,
                                                ref_buffer_,
                                                max_errors,
                                                scratch_banded_prev_,
                                                scratch_banded_curr_);
            best_score = std::min(best_score, score);
        }
        if (best_score == 0) {
            return 0;
        }
    }
    return best_score;
}

AlignmentHit MapperEngine::verify_candidate(const std::string& oriented_read,
                                            const CandidateInfo& candidate,
                                            int max_errors) {
    AlignmentHit no_hit;
    no_hit.edit_distance = std::numeric_limits<int>::max();

    const uint32_t global_start = candidate.start;
    const std::size_t chrom_index = index_.chromosome_for_position(global_start);
    const auto& chrom = index_.chromosome(chrom_index);
    const uint64_t local_pos = static_cast<uint64_t>(global_start) - chrom.start;
    const uint32_t min_ref_len = oriented_read.size() > static_cast<std::size_t>(max_errors)
        ? static_cast<uint32_t>(oriented_read.size() - static_cast<std::size_t>(max_errors))
        : 1u;
    const uint32_t max_ref_len = std::min<uint32_t>(
        static_cast<uint32_t>(oriented_read.size() + static_cast<std::size_t>(max_errors)),
        chrom.length - static_cast<uint32_t>(local_pos));

    if (min_ref_len > max_ref_len) {
        return no_hit;
    }

    AlignmentHit best_hit = no_hit;
    for (uint32_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        index_.extract_sequence(global_start, ref_len, ref_buffer_);
        AlignmentResult aligned = alignment_workspace_.align(oriented_read, ref_buffer_, max_errors);
        if (aligned.edit_distance > max_errors) {
            continue;
        }

        AlignmentHit candidate;
        candidate.chrom_index = chrom_index;
        candidate.global_pos = global_start;
        candidate.ref_pos_1based = static_cast<uint32_t>(local_pos + 1u);
        candidate.edit_distance = aligned.edit_distance;
        candidate.cigar = std::move(aligned.cigar);
        if (best_hit.edit_distance == std::numeric_limits<int>::max() || better_hit(candidate, best_hit)) {
            best_hit = std::move(candidate);
        }
    }

    return best_hit;
}

bool MapperEngine::better_hit(const AlignmentHit& lhs, const AlignmentHit& rhs) const {
    if (lhs.edit_distance != rhs.edit_distance) {
        return lhs.edit_distance < rhs.edit_distance;
    }
    const auto& lchrom = index_.chromosome(lhs.chrom_index).name;
    const auto& rchrom = index_.chromosome(rhs.chrom_index).name;
    if (lchrom != rchrom) {
        return lchrom < rchrom;
    }
    if (lhs.ref_pos_1based != rhs.ref_pos_1based) {
        return lhs.ref_pos_1based < rhs.ref_pos_1based;
    }
    return lhs.cigar < rhs.cigar;
}

void MapperEngine::maybe_add_hit(std::vector<AlignmentHit>& hits, AlignmentHit&& hit) const {
    if (hit.edit_distance == std::numeric_limits<int>::max()) {
        return;
    }
    for (AlignmentHit& existing : hits) {
        if (existing.chrom_index == hit.chrom_index && existing.global_pos == hit.global_pos && existing.cigar == hit.cigar) {
            if (better_hit(hit, existing)) {
                existing = std::move(hit);
            }
            return;
        }
    }
    hits.push_back(std::move(hit));
    std::sort(hits.begin(), hits.end(), [this](const AlignmentHit& lhs, const AlignmentHit& rhs) {
        return better_hit(lhs, rhs);
    });
    if (hits.size() > 2u) {
        hits.resize(2u);
    }
}

void MapperEngine::search_orientation(const std::string& oriented_read,
                                      int max_errors,
                                      std::vector<AlignmentHit>& hits) {
    scratch_seeds_.clear();
    scratch_prefiltered_.clear();
    collect_seeds(oriented_read, max_errors, scratch_seeds_);
    if (scratch_seeds_.empty()) {
        return;
    }

    const MyersQuery query = build_myers_query(oriented_read);
    scratch_candidates_.clear();
    generate_candidates(scratch_seeds_, oriented_read.size(), max_errors, scratch_candidates_);

    for (const CandidateInfo& candidate : scratch_candidates_) {
        const int best_score = prefilter_candidate(oriented_read, query, candidate, max_errors);
        if (best_score <= max_errors) {
            scratch_prefiltered_.push_back(PrefilterCandidate{candidate, best_score});
        }
    }

    std::sort(scratch_prefiltered_.begin(), scratch_prefiltered_.end(),
              [](const PrefilterCandidate& lhs, const PrefilterCandidate& rhs) {
                  if (lhs.best_score != rhs.best_score) {
                      return lhs.best_score < rhs.best_score;
                  }
                  if (lhs.candidate.support != rhs.candidate.support) {
                      return lhs.candidate.support > rhs.candidate.support;
                  }
                  if (lhs.candidate.best_seed_freq != rhs.candidate.best_seed_freq) {
                      return lhs.candidate.best_seed_freq < rhs.candidate.best_seed_freq;
                  }
                  return lhs.candidate.start < rhs.candidate.start;
              });

    const std::size_t dp_limit = std::min<std::size_t>(kMaxDpFinalistsPerOrientation, scratch_prefiltered_.size());
    for (std::size_t i = 0; i < dp_limit; ++i) {
        maybe_add_hit(hits, verify_candidate(oriented_read, scratch_prefiltered_[i].candidate, max_errors));
    }
}

std::string MapperEngine::format_record(const FastqRecord& record,
                                        const std::vector<AlignmentHit>& hits) const {
    std::string line;
    line.reserve(256);
    if (hits.empty()) {
        line += record.name;
        line += " * 0 * ";
        line += record.seq;
        line += ' ';
        line += record.qual;
        line += '\n';
        return line;
    }

    const AlignmentHit& primary = hits.front();
    line += record.name;
    line += ' ';
    line += index_.chromosome(primary.chrom_index).name;
    line += ' ';
    line += std::to_string(primary.ref_pos_1based);
    line += ' ';
    line += primary.cigar;
    line += ' ';
    line += record.seq;
    line += ' ';
    line += record.qual;
    if (hits.size() > 1u) {
        const AlignmentHit& alt = hits[1];
        line += " ALT:";
        line += index_.chromosome(alt.chrom_index).name;
        line += ',';
        line += std::to_string(alt.ref_pos_1based);
        line += ',';
        line += alt.cigar;
    }
    line += '\n';
    return line;
}

std::string MapperEngine::map_record(const FastqRecord& record, int max_errors) {
    const std::string normalized = normalize_sequence(record.seq);
    const std::string revcomp = reverse_complement(normalized);
    std::vector<AlignmentHit> hits;
    hits.reserve(2);

    search_orientation(normalized, max_errors, hits);
    search_orientation(revcomp, max_errors, hits);
    std::sort(hits.begin(), hits.end(), [this](const AlignmentHit& lhs, const AlignmentHit& rhs) {
        return better_hit(lhs, rhs);
    });
    if (hits.size() > 2u) {
        hits.resize(2u);
    }
    return format_record(record, hits);
}

}  // namespace mapper_speed

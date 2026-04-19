#include "search.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <string_view>

#include "common.hpp"
#include "nucleotide.hpp"

namespace mapper_speed {

namespace {

constexpr std::size_t kMaxDpFinalistsPerOrientation = 10;
constexpr std::size_t kMaxCandidatePrefilterPerOrientation = 20;
constexpr uint32_t kSeedMinSpacing = 10;
constexpr uint32_t kNearDuplicateTolerance = 10;
constexpr uint32_t kCandidateClusterSlack = 3;
constexpr uint32_t kQuickSeedFreq = 8;
constexpr std::size_t kQuickSeedCount = 3;
constexpr std::size_t kQuickCandidateLimit = 8;

std::size_t target_seed_count(int max_errors) {
    switch (max_errors) {
        case 0: return 10;
        case 1: return 12;
        case 2: return 14;
        default: return 16;
    }
}

bool far_enough_from_selected(const std::vector<SeedSpec>& selected,
                              uint32_t read_offset,
                              uint32_t min_spacing) {
    for (const SeedSpec& seed : selected) {
        const uint32_t delta = (seed.read_offset > read_offset)
            ? (seed.read_offset - read_offset)
            : (read_offset - seed.read_offset);
        if (delta < min_spacing) {
            return false;
        }
    }
    return true;
}

}  // namespace

MapperEngine::MapperEngine(const IndexView& index)
    : index_(index),
      dispatch_(resolve_myers_dispatch()) {
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
    const uint32_t stride = std::max<uint32_t>(1u, index_.index_stride());
    if (target <= 1u) {
        positions.push_back(0);
        return;
    }

    positions.reserve(target + stride);
    const uint32_t residue_limit = std::min<uint32_t>(stride, max_start + 1u);
    for (uint32_t residue = 0; residue < residue_limit; ++residue) {
        positions.push_back(residue);
    }
    for (std::size_t i = 0; i < target; ++i) {
        const uint64_t numerator = static_cast<uint64_t>(i) * max_start + (target - 1u) / 2u;
        positions.push_back(static_cast<uint32_t>(numerator / (target - 1u)));
    }
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
}

void MapperEngine::collect_seeds(const std::string& normalized_read,
                                 int max_errors,
                                 std::vector<SeedSpec>& seeds) {
    seeds.clear();
    scratch_seed_candidates_.clear();
    scratch_seed_positions_.clear();
    collect_seed_positions(normalized_read.size(), max_errors, scratch_seed_positions_);
    if (scratch_seed_positions_.empty()) {
        return;
    }

    const uint32_t max_start = static_cast<uint32_t>(normalized_read.size() - kSeedLength);
    scratch_seed_keys_.assign(max_start + 1u, 0u);
    scratch_seed_valid_.assign(max_start + 1u, 0u);
    uint32_t rolling_key = 0;
    uint32_t valid_run = 0;
    for (uint32_t i = 0; i < normalized_read.size(); ++i) {
        const uint8_t code = base_to_code(normalized_read[i]);
        if (code == kBaseCodeInvalid) {
            rolling_key = 0;
            valid_run = 0;
            continue;
        }
        rolling_key = (rolling_key << 2u) | code;
        if (valid_run < kSeedLength) {
            ++valid_run;
        }
        if (valid_run >= kSeedLength) {
            const uint32_t start = i - kSeedLength + 1u;
            scratch_seed_keys_[start] = rolling_key;
            scratch_seed_valid_[start] = 1u;
        }
    }

    for (uint32_t pos : scratch_seed_positions_) {
        if (pos > max_start || scratch_seed_valid_[pos] == 0u) {
            continue;
        }
        const uint32_t key = scratch_seed_keys_[pos];
        const uint32_t freq = index_.occurrence_count(key);
        if (freq == 0 || freq > kHighFreqAllowFallback) {
            continue;
        }
        scratch_seed_candidates_.push_back(SeedSpec{pos, key, freq});
    }

    if (scratch_seed_candidates_.empty()) {
        return;
    }

    std::sort(scratch_seed_candidates_.begin(), scratch_seed_candidates_.end(),
              [](const SeedSpec& lhs, const SeedSpec& rhs) {
                  if (lhs.frequency != rhs.frequency) {
                      return lhs.frequency < rhs.frequency;
                  }
                  return lhs.read_offset < rhs.read_offset;
              });

    const std::size_t target = target_seed_count(max_errors);
    for (uint32_t spacing : {kSeedMinSpacing, kSeedMinSpacing / 2u, 0u}) {
        for (const SeedSpec& candidate : scratch_seed_candidates_) {
            if (candidate.frequency > kHighFreqSkip) {
                continue;
            }
            if (!far_enough_from_selected(seeds, candidate.read_offset, spacing)) {
                continue;
            }
            seeds.push_back(candidate);
            if (seeds.size() >= target) {
                break;
            }
        }
        if (seeds.size() >= target) {
            break;
        }
    }

    for (uint32_t spacing : {kSeedMinSpacing, kSeedMinSpacing / 2u, 0u}) {
        for (const SeedSpec& candidate : scratch_seed_candidates_) {
            if (!far_enough_from_selected(seeds, candidate.read_offset, spacing)) {
                continue;
            }
            seeds.push_back(candidate);
            if (seeds.size() >= target) {
                break;
            }
        }
        if (seeds.size() >= target) {
            break;
        }
    }

    std::sort(seeds.begin(), seeds.end(), [](const SeedSpec& lhs, const SeedSpec& rhs) {
        if (lhs.frequency != rhs.frequency) {
            return lhs.frequency < rhs.frequency;
        }
        return lhs.read_offset < rhs.read_offset;
    });
    seeds.erase(std::unique(seeds.begin(), seeds.end(),
                            [](const SeedSpec& lhs, const SeedSpec& rhs) {
                                return lhs.read_offset == rhs.read_offset && lhs.key == rhs.key;
                            }),
                seeds.end());
}

bool MapperEngine::try_upfront_exactish(const std::string& oriented_read,
                                        int max_errors,
                                        const std::vector<SeedSpec>& seeds,
                                        std::vector<AlignmentHit>& hits) {
    if (max_errors > 1 || seeds.empty()) {
        return false;
    }

    scratch_candidates_.clear();
    const std::size_t seed_limit = std::min<std::size_t>(kQuickSeedCount, seeds.size());
    for (std::size_t i = 0; i < seed_limit; ++i) {
        const SeedSpec& seed = seeds[i];
        if (seed.frequency > kQuickSeedFreq) {
            break;
        }
        const auto range = index_.positions_for(seed.key);
        for (const uint32_t* it = range.first; it != range.second; ++it) {
            const int64_t candidate_start = static_cast<int64_t>(*it) - static_cast<int64_t>(seed.read_offset);
            if (candidate_start < 0) {
                continue;
            }
            const uint32_t start = static_cast<uint32_t>(candidate_start);
            const uint16_t chrom_index = static_cast<uint16_t>(index_.chromosome_for_position(start));
            if (!index_.stays_within_chromosome(start, static_cast<uint32_t>(oriented_read.size()))) {
                continue;
            }
            scratch_candidates_.push_back(CandidateInfo{start, seed.frequency, 1, chrom_index});
        }
    }

    if (scratch_candidates_.empty()) {
        return false;
    }

    std::sort(scratch_candidates_.begin(), scratch_candidates_.end(),
              [](const CandidateInfo& lhs, const CandidateInfo& rhs) {
                  if (lhs.start != rhs.start) {
                      return lhs.start < rhs.start;
                  }
                  return lhs.best_seed_freq < rhs.best_seed_freq;
              });
    scratch_candidates_.erase(std::unique(scratch_candidates_.begin(), scratch_candidates_.end(),
                                          [](const CandidateInfo& lhs, const CandidateInfo& rhs) {
                                              return lhs.start == rhs.start;
                                          }),
                              scratch_candidates_.end());

    const std::size_t limit = std::min<std::size_t>(kQuickCandidateLimit, scratch_candidates_.size());
    for (std::size_t i = 0; i < limit; ++i) {
        const CandidateInfo& candidate = scratch_candidates_[i];
        const std::size_t chrom_index = candidate.chrom_index;
        const auto& chrom = index_.chromosome(chrom_index);
        const uint64_t local_pos = static_cast<uint64_t>(candidate.start) - chrom.start;
        index_.extract_sequence(candidate.start, static_cast<uint32_t>(oriented_read.size()), ref_buffer_);
        const uint32_t mismatches =
            count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), oriented_read.size());
        if (static_cast<int>(mismatches) > max_errors) {
            continue;
        }

        AlignmentHit hit;
        hit.chrom_index = chrom_index;
        hit.global_pos = candidate.start;
        hit.ref_pos_1based = static_cast<uint32_t>(local_pos + 1u);
        hit.edit_distance = static_cast<int>(mismatches);
        hit.cigar = std::to_string(oriented_read.size()) + "M";
        maybe_add_hit(hits, std::move(hit));

        if (mismatches == 0u) {
            return true;
        }
        if (hits.size() >= 2u && hits.back().edit_distance <= 1) {
            return true;
        }
    }

    return false;
}

void MapperEngine::generate_candidates(const std::vector<SeedSpec>& seeds,
                                       std::size_t read_length,
                                       int max_errors,
                                       std::vector<CandidateInfo>& starts) {
    starts.clear();
    scratch_anchors_.clear();
    const uint32_t cluster_window = static_cast<uint32_t>(max_errors) + kCandidateClusterSlack;
    const uint32_t min_ref_len = read_length > static_cast<std::size_t>(max_errors)
        ? static_cast<uint32_t>(read_length - static_cast<std::size_t>(max_errors))
        : 1u;

    for (uint16_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        const SeedSpec& seed = seeds[seed_index];
        const auto range = index_.positions_for(seed.key);
        for (const uint32_t* it = range.first; it != range.second; ++it) {
            const int64_t base_candidate = static_cast<int64_t>(*it) - static_cast<int64_t>(seed.read_offset);
            for (int delta = -max_errors; delta <= max_errors; ++delta) {
                const int64_t candidate = base_candidate + delta;
                if (candidate < 0) {
                    continue;
                }
                const uint32_t start = static_cast<uint32_t>(candidate);
                if (static_cast<uint64_t>(start) + min_ref_len > index_.genome_length()) {
                    continue;
                }
                scratch_anchors_.push_back(SeedAnchor{start, seed_index, 0});
            }
        }
    }

    if (scratch_anchors_.empty()) {
        return;
    }

    std::sort(scratch_anchors_.begin(), scratch_anchors_.end(),
              [](const SeedAnchor& lhs, const SeedAnchor& rhs) {
                  return lhs.start < rhs.start;
              });

    struct ChainCluster {
        uint32_t min_start = 0;
        uint32_t max_start = 0;
        uint32_t best_seed_freq = 0xFFFFFFFFu;
        uint32_t seed_mask = 0;
        uint16_t anchor_count = 0;
        uint16_t unique_seed_count = 0;
        uint32_t min_read_offset = 0xFFFFFFFFu;
        uint32_t max_read_offset = 0;
        uint16_t chrom_index = 0;
    };

    std::vector<ChainCluster> clusters;
    clusters.reserve(scratch_anchors_.size());

    for (const SeedAnchor& anchor : scratch_anchors_) {
        const uint16_t chrom_index = static_cast<uint16_t>(index_.chromosome_for_position(anchor.start));
        const SeedSpec& seed = seeds[anchor.seed_index];
        if (!clusters.empty()) {
            ChainCluster& cluster = clusters.back();
            if (cluster.chrom_index == chrom_index &&
                anchor.start >= cluster.max_start &&
                anchor.start - cluster.max_start <= cluster_window) {
                cluster.max_start = anchor.start;
                cluster.best_seed_freq = std::min(cluster.best_seed_freq, seed.frequency);
                cluster.anchor_count = static_cast<uint16_t>(std::min<uint32_t>(cluster.anchor_count + 1u, 0xFFFFu));
                cluster.min_read_offset = std::min(cluster.min_read_offset, seed.read_offset);
                cluster.max_read_offset = std::max(cluster.max_read_offset, seed.read_offset);
                const uint32_t bit = (anchor.seed_index < 32u) ? (uint32_t{1} << anchor.seed_index) : 0u;
                if ((cluster.seed_mask & bit) == 0u) {
                    cluster.seed_mask |= bit;
                    ++cluster.unique_seed_count;
                }
                continue;
            }
        }

        ChainCluster cluster{};
        cluster.min_start = anchor.start;
        cluster.max_start = anchor.start;
        cluster.best_seed_freq = seed.frequency;
        cluster.anchor_count = 1;
        cluster.unique_seed_count = 1;
        cluster.min_read_offset = seed.read_offset;
        cluster.max_read_offset = seed.read_offset;
        cluster.chrom_index = chrom_index;
        cluster.seed_mask = (anchor.seed_index < 32u) ? (uint32_t{1} << anchor.seed_index) : 0u;
        clusters.push_back(cluster);
    }

    std::sort(clusters.begin(), clusters.end(), [](const ChainCluster& lhs, const ChainCluster& rhs) {
        if (lhs.unique_seed_count != rhs.unique_seed_count) {
            return lhs.unique_seed_count > rhs.unique_seed_count;
        }
        const uint32_t lhs_span = lhs.max_read_offset - lhs.min_read_offset;
        const uint32_t rhs_span = rhs.max_read_offset - rhs.min_read_offset;
        if (lhs_span != rhs_span) {
            return lhs_span > rhs_span;
        }
        if (lhs.anchor_count != rhs.anchor_count) {
            return lhs.anchor_count > rhs.anchor_count;
        }
        if (lhs.best_seed_freq != rhs.best_seed_freq) {
            return lhs.best_seed_freq < rhs.best_seed_freq;
        }
        return lhs.min_start < rhs.min_start;
    });

    const std::size_t limit = std::min<std::size_t>(kMaxVerifyPerOrientation, clusters.size());
    for (std::size_t i = 0; i < limit; ++i) {
        const uint32_t representative = clusters[i].min_start + ((clusters[i].max_start - clusters[i].min_start) / 2u);
        starts.push_back(CandidateInfo{
            representative,
            clusters[i].best_seed_freq,
            clusters[i].unique_seed_count,
            static_cast<uint16_t>(clusters[i].chrom_index)
        });
    }
}

int MapperEngine::prefilter_candidate(const std::string& oriented_read,
                                      const MyersQuery& query,
                                      const CandidateInfo& candidate,
                                      int max_errors) {
    const uint32_t global_start = candidate.start;
    const std::size_t chrom_index = candidate.chrom_index;
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
            if (mismatches <= 1u) {
                return static_cast<int>(mismatches);
            }
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
    const std::size_t chrom_index = candidate.chrom_index;
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
    if (index_.stays_within_chromosome(global_start, static_cast<uint32_t>(oriented_read.size()))) {
        index_.extract_sequence(global_start, static_cast<uint32_t>(oriented_read.size()), ref_buffer_);
        const uint32_t mismatches =
            count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), oriented_read.size());
        if (mismatches <= 1u) {
            best_hit.chrom_index = chrom_index;
            best_hit.global_pos = global_start;
            best_hit.ref_pos_1based = static_cast<uint32_t>(local_pos + 1u);
            best_hit.edit_distance = static_cast<int>(mismatches);
            best_hit.cigar = std::to_string(oriented_read.size()) + "M";
            return best_hit;
        }
    }
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
        if (existing.chrom_index == hit.chrom_index && existing.edit_distance == hit.edit_distance) {
            const uint32_t delta = (existing.ref_pos_1based > hit.ref_pos_1based)
                ? (existing.ref_pos_1based - hit.ref_pos_1based)
                : (hit.ref_pos_1based - existing.ref_pos_1based);
            if (delta <= kNearDuplicateTolerance) {
                if (better_hit(hit, existing)) {
                    existing = std::move(hit);
                }
                return;
            }
        }
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
    if (try_upfront_exactish(oriented_read, max_errors, scratch_seeds_, hits)) {
        return;
    }

    const MyersQuery query = build_myers_query(oriented_read);
    scratch_candidates_.clear();
    generate_candidates(scratch_seeds_, oriented_read.size(), max_errors, scratch_candidates_);

    const int initial_score_limit = (hits.size() >= 2u) ? hits.back().edit_distance : max_errors;
    for (const CandidateInfo& candidate : scratch_candidates_) {
        const int best_score = prefilter_candidate(oriented_read, query, candidate, max_errors);
        if (best_score <= initial_score_limit) {
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

    const int best_prefilter_score = scratch_prefiltered_.empty()
        ? std::numeric_limits<int>::max()
        : scratch_prefiltered_.front().best_score;
    const std::size_t dp_limit =
        std::min<std::size_t>(kMaxDpFinalistsPerOrientation,
                              std::min<std::size_t>(kMaxCandidatePrefilterPerOrientation, scratch_prefiltered_.size()));
    for (std::size_t i = 0; i < dp_limit; ++i) {
        const int current_score_limit = (hits.size() >= 2u) ? hits.back().edit_distance : max_errors;
        if (scratch_prefiltered_[i].best_score > current_score_limit ||
            scratch_prefiltered_[i].best_score > best_prefilter_score + 1) {
            break;
        }
        maybe_add_hit(hits, verify_candidate(oriented_read, scratch_prefiltered_[i].candidate, max_errors));
        if (hits.size() >= 2u && hits.back().edit_distance == 0) {
            break;
        }
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
    std::string normalized;
    std::string revcomp;
    normalize_and_reverse_complement(record.seq, normalized, revcomp);
    std::vector<AlignmentHit> hits;
    hits.reserve(2);

    search_orientation(normalized, max_errors, hits);
    if (hits.size() < 2u || hits.back().edit_distance > 0) {
        search_orientation(revcomp, max_errors, hits);
    }
    std::sort(hits.begin(), hits.end(), [this](const AlignmentHit& lhs, const AlignmentHit& rhs) {
        return better_hit(lhs, rhs);
    });
    if (hits.size() > 2u) {
        hits.resize(2u);
    }
    return format_record(record, hits);
}

}  // namespace mapper_speed

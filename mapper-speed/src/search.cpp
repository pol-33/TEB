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
constexpr std::size_t kDenseFullSeedProbeLimit = 192;
constexpr std::size_t kDenseSampledSeedProbeCount = 48;
constexpr std::size_t kSeedOccurrenceExpandLimitPrimary = 1024;
constexpr std::size_t kSeedOccurrenceExpandLimitSecondary = 256;
constexpr uint32_t kSeedMinSpacing = 10;
constexpr uint32_t kNearDuplicateTolerance = 10;
constexpr uint32_t kCandidateClusterSlack = 3;
constexpr uint32_t kQuickSeedFreq = 8;
constexpr std::size_t kQuickSeedCount = 3;
constexpr std::size_t kQuickCandidateLimit = 8;
constexpr std::size_t kCandidatePrefetchDistance = 4;

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

    if (stride == 1u) {
        const std::size_t position_count = static_cast<std::size_t>(max_start) + 1u;
        if (position_count <= kDenseFullSeedProbeLimit) {
            positions.reserve(position_count);
            for (uint32_t pos = 0; pos <= max_start; ++pos) {
                positions.push_back(pos);
            }
            return;
        }

        const std::size_t probe_target =
            std::min<std::size_t>(position_count, std::max<std::size_t>(target * 4u, kDenseSampledSeedProbeCount));
        positions.reserve(probe_target);
        for (std::size_t i = 0; i < probe_target; ++i) {
            const uint64_t numerator = static_cast<uint64_t>(i) * max_start + (probe_target - 1u) / 2u;
            positions.push_back(static_cast<uint32_t>(numerator / (probe_target - 1u)));
        }
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
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
        const auto range = index_.positions_for(key);
        const uint32_t freq = static_cast<uint32_t>(range.second - range.first);
        if (freq == 0 || freq > kHighFreqAllowFallback) {
            continue;
        }
        scratch_seed_candidates_.push_back(SeedSpec{pos, key, freq, range.first, range.second});
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
                                        const std::vector<uint8_t>& packed_read,
                                        bool packed_read_valid,
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
        for (const uint32_t* it = seed.begin; it != seed.end; ++it) {
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
        if (i + kCandidatePrefetchDistance < limit) {
            index_.prefetch_reference_window(scratch_candidates_[i + kCandidatePrefetchDistance].start,
                                             static_cast<uint32_t>(oriented_read.size() + max_errors));
        }
        const CandidateInfo& candidate = scratch_candidates_[i];
        const std::size_t chrom_index = candidate.chrom_index;
        const auto& chrom = index_.chromosome(chrom_index);
        const uint64_t local_pos = static_cast<uint64_t>(candidate.start) - chrom.start;
        uint32_t mismatches = 0;
        if (packed_read_valid) {
            mismatches = index_.count_mismatches_packed(candidate.start,
                                                       packed_read.data(),
                                                       static_cast<uint32_t>(oriented_read.size()),
                                                       static_cast<uint32_t>(max_errors));
        } else {
            index_.extract_sequence(candidate.start, static_cast<uint32_t>(oriented_read.size()), ref_buffer_);
            mismatches =
                count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), oriented_read.size());
        }
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
    scratch_clusters_.clear();
    const uint32_t cluster_window = static_cast<uint32_t>(max_errors) + kCandidateClusterSlack;
    const uint32_t min_ref_len = read_length > static_cast<std::size_t>(max_errors)
        ? static_cast<uint32_t>(read_length - static_cast<std::size_t>(max_errors))
        : 1u;
    const bool dense_index = index_.index_stride() == 1u;
    const std::size_t delta_count = static_cast<std::size_t>((max_errors * 2) + 1);

    auto occurrence_expand_limit = [&](std::size_t seed_rank, uint32_t frequency) -> std::size_t {
        std::size_t limit = (seed_rank < 2u) ? kSeedOccurrenceExpandLimitPrimary : kSeedOccurrenceExpandLimitSecondary;
        if (dense_index) {
            limit = std::max<std::size_t>(64u, limit / 2u);
        }
        if (frequency <= 32u) {
            limit = std::max<std::size_t>(limit, static_cast<std::size_t>(frequency));
        }
        return limit;
    };

    std::size_t estimated_anchor_count = 0;
    for (std::size_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        const SeedSpec& seed = seeds[seed_index];
        const std::size_t occurrence_count = static_cast<std::size_t>(seed.end - seed.begin);
        const std::size_t expand_limit = occurrence_expand_limit(seed_index, seed.frequency);
        estimated_anchor_count += std::min(occurrence_count, expand_limit) * delta_count;
    }
    scratch_anchors_.reserve(std::max(scratch_anchors_.capacity(), estimated_anchor_count));

    auto expand_occurrence = [&](uint16_t seed_index, uint16_t chrom_index, uint32_t hit_pos) {
        const SeedSpec& seed = seeds[seed_index];
        const auto& chrom = index_.chromosome(chrom_index);
        const int64_t base_candidate = static_cast<int64_t>(hit_pos) - static_cast<int64_t>(seed.read_offset);
        const int64_t min_start = static_cast<int64_t>(chrom.start);
        const int64_t max_start =
            static_cast<int64_t>(chrom.start) + static_cast<int64_t>(chrom.length) - static_cast<int64_t>(min_ref_len);
        const int64_t delta_min = std::max<int64_t>(-max_errors, min_start - base_candidate);
        const int64_t delta_max = std::min<int64_t>(max_errors, max_start - base_candidate);
        if (delta_min > delta_max) {
            return;
        }
        const uint32_t seed_bit = (seed_index < 32u) ? (uint32_t{1} << seed_index) : 0u;
        for (int64_t delta = delta_min; delta <= delta_max; ++delta) {
            scratch_anchors_.push_back(SeedAnchor{
                static_cast<uint32_t>(base_candidate + delta),
                seed.frequency,
                seed_bit,
                seed.read_offset,
                seed.read_offset,
                1,
                1,
                chrom_index
            });
        }
    };

    for (uint16_t seed_index = 0; seed_index < seeds.size(); ++seed_index) {
        const SeedSpec& seed = seeds[seed_index];
        const std::size_t occurrence_count = static_cast<std::size_t>(seed.end - seed.begin);
        const std::size_t expand_limit = occurrence_expand_limit(seed_index, seed.frequency);
        if (occurrence_count <= expand_limit) {
            for (const uint32_t* it = seed.begin; it != seed.end; ++it) {
                if (static_cast<std::size_t>(seed.end - it) > 8u) {
                    prefetch_read_mostly(it + 8);
                }
                const uint16_t chrom_index =
                    static_cast<uint16_t>(index_.chromosome_for_position(*it));
                expand_occurrence(seed_index, chrom_index, *it);
            }
            continue;
        }

        for (std::size_t i = 0; i < expand_limit; ++i) {
            const uint64_t numerator =
                static_cast<uint64_t>(i) * static_cast<uint64_t>(occurrence_count - 1u) +
                static_cast<uint64_t>(expand_limit - 1u) / 2u;
            const std::size_t sample_index =
                static_cast<std::size_t>(numerator / static_cast<uint64_t>(expand_limit - 1u));
            if (i + 2u < expand_limit) {
                const uint64_t next_numerator =
                    static_cast<uint64_t>(i + 2u) * static_cast<uint64_t>(occurrence_count - 1u) +
                    static_cast<uint64_t>(expand_limit - 1u) / 2u;
                const std::size_t next_sample_index =
                    static_cast<std::size_t>(next_numerator / static_cast<uint64_t>(expand_limit - 1u));
                prefetch_read_mostly(seed.begin + next_sample_index);
            }
            const uint32_t hit_pos = seed.begin[sample_index];
            const uint16_t chrom_index =
                static_cast<uint16_t>(index_.chromosome_for_position(hit_pos));
            expand_occurrence(seed_index, chrom_index, hit_pos);
        }
    }

    if (scratch_anchors_.empty()) {
        return;
    }

    std::sort(scratch_anchors_.begin(), scratch_anchors_.end(),
              [](const SeedAnchor& lhs, const SeedAnchor& rhs) {
                  if (lhs.chrom_index != rhs.chrom_index) {
                      return lhs.chrom_index < rhs.chrom_index;
                  }
                  if (lhs.start != rhs.start) {
                      return lhs.start < rhs.start;
                  }
                  return lhs.best_seed_freq < rhs.best_seed_freq;
              });

    std::size_t write_index = 0;
    for (std::size_t i = 0; i < scratch_anchors_.size(); ++i) {
        if (write_index == 0u ||
            scratch_anchors_[i].chrom_index != scratch_anchors_[write_index - 1u].chrom_index ||
            scratch_anchors_[i].start != scratch_anchors_[write_index - 1u].start) {
            if (write_index != i) {
                scratch_anchors_[write_index] = scratch_anchors_[i];
            }
            ++write_index;
            continue;
        }

        SeedAnchor& merged = scratch_anchors_[write_index - 1u];
        const SeedAnchor& current = scratch_anchors_[i];
        merged.best_seed_freq = std::min(merged.best_seed_freq, current.best_seed_freq);
        merged.anchor_count = static_cast<uint16_t>(std::min<uint32_t>(merged.anchor_count + current.anchor_count, 0xFFFFu));
        merged.min_read_offset = std::min(merged.min_read_offset, current.min_read_offset);
        merged.max_read_offset = std::max(merged.max_read_offset, current.max_read_offset);
        const uint32_t new_bits = current.seed_mask & ~merged.seed_mask;
        merged.seed_mask |= current.seed_mask;
        merged.unique_seed_count = static_cast<uint16_t>(
            std::min<uint32_t>(merged.unique_seed_count + static_cast<uint32_t>(__builtin_popcount(new_bits)),
                               0xFFFFu));
    }
    scratch_anchors_.resize(write_index);

    scratch_clusters_.reserve(std::max(scratch_clusters_.capacity(), scratch_anchors_.size()));

    for (const SeedAnchor& anchor : scratch_anchors_) {
        const uint16_t chrom_index = anchor.chrom_index;
        if (!scratch_clusters_.empty()) {
            CandidateCluster& cluster = scratch_clusters_.back();
            if (cluster.chrom_index == chrom_index &&
                anchor.start >= cluster.max_start &&
                anchor.start - cluster.max_start <= cluster_window) {
                cluster.max_start = anchor.start;
                cluster.best_seed_freq = std::min(cluster.best_seed_freq, anchor.best_seed_freq);
                cluster.anchor_count = static_cast<uint16_t>(
                    std::min<uint32_t>(cluster.anchor_count + anchor.anchor_count, 0xFFFFu));
                cluster.min_read_offset = std::min(cluster.min_read_offset, anchor.min_read_offset);
                cluster.max_read_offset = std::max(cluster.max_read_offset, anchor.max_read_offset);
                const uint32_t new_bits = anchor.seed_mask & ~cluster.seed_mask;
                cluster.seed_mask |= anchor.seed_mask;
                cluster.unique_seed_count = static_cast<uint16_t>(
                    std::min<uint32_t>(cluster.unique_seed_count + static_cast<uint32_t>(__builtin_popcount(new_bits)),
                                       0xFFFFu));
                continue;
            }
        }

        CandidateCluster cluster{};
        cluster.min_start = anchor.start;
        cluster.max_start = anchor.start;
        cluster.best_seed_freq = anchor.best_seed_freq;
        cluster.anchor_count = anchor.anchor_count;
        cluster.unique_seed_count = anchor.unique_seed_count;
        cluster.min_read_offset = anchor.min_read_offset;
        cluster.max_read_offset = anchor.max_read_offset;
        cluster.chrom_index = chrom_index;
        cluster.seed_mask = anchor.seed_mask;
        scratch_clusters_.push_back(cluster);
    }

    const auto cluster_better = [](const CandidateCluster& lhs, const CandidateCluster& rhs) {
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
    };

    const std::size_t limit = std::min<std::size_t>(kMaxVerifyPerOrientation, scratch_clusters_.size());
    if (limit == 0u) {
        return;
    }
    starts.reserve(limit);
    std::partial_sort(scratch_clusters_.begin(), scratch_clusters_.begin() + limit, scratch_clusters_.end(),
                      cluster_better);
    for (std::size_t i = 0; i < limit; ++i) {
        const uint32_t representative =
            scratch_clusters_[i].min_start + ((scratch_clusters_[i].max_start - scratch_clusters_[i].min_start) / 2u);
        starts.push_back(CandidateInfo{
            representative,
            scratch_clusters_[i].best_seed_freq,
            scratch_clusters_[i].unique_seed_count,
            static_cast<uint16_t>(scratch_clusters_[i].chrom_index)
        });
    }
}

int MapperEngine::prefilter_candidate(const std::string& oriented_read,
                                      const std::vector<uint8_t>& packed_read,
                                      bool packed_read_valid,
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
        if (ref_len == oriented_read.size()) {
            uint32_t mismatches = 0;
            if (packed_read_valid) {
                mismatches = index_.count_mismatches_packed(global_start,
                                                           packed_read.data(),
                                                           ref_len,
                                                           static_cast<uint32_t>(max_errors));
            } else {
                index_.extract_sequence(global_start, ref_len, ref_buffer_);
                mismatches = count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), ref_len);
            }
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
            if (packed_read_valid) {
                index_.extract_sequence(global_start, ref_len, ref_buffer_);
            }
            const int score = bounded_edit_distance(dispatch_, query, ref_buffer_, max_errors);
            best_score = std::min(best_score, score);
        } else {
            index_.extract_sequence(global_start, ref_len, ref_buffer_);
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
                                            const std::vector<uint8_t>& packed_read,
                                            bool packed_read_valid,
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
        uint32_t mismatches = 0;
        if (packed_read_valid) {
            mismatches = index_.count_mismatches_packed(global_start,
                                                       packed_read.data(),
                                                       static_cast<uint32_t>(oriented_read.size()),
                                                       1u);
        } else {
            index_.extract_sequence(global_start, static_cast<uint32_t>(oriented_read.size()), ref_buffer_);
            mismatches =
                count_byte_mismatches_bulk(oriented_read.data(), ref_buffer_.data(), oriented_read.size());
        }
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
                                      const std::vector<uint8_t>& packed_read,
                                      bool packed_read_valid,
                                      int max_errors,
                                      std::vector<AlignmentHit>& hits) {
    scratch_seeds_.clear();
    scratch_prefiltered_.clear();
    collect_seeds(oriented_read, max_errors, scratch_seeds_);
    if (scratch_seeds_.empty()) {
        return;
    }
    if (try_upfront_exactish(oriented_read, packed_read, packed_read_valid, max_errors, scratch_seeds_, hits)) {
        return;
    }

    const MyersQuery query = build_myers_query(oriented_read);
    scratch_candidates_.clear();
    generate_candidates(scratch_seeds_, oriented_read.size(), max_errors, scratch_candidates_);

    const int initial_score_limit = (hits.size() >= 2u) ? hits.back().edit_distance : max_errors;
    scratch_prefiltered_.reserve(scratch_candidates_.size());
    for (const CandidateInfo& candidate : scratch_candidates_) {
        const std::size_t candidate_index = static_cast<std::size_t>(&candidate - scratch_candidates_.data());
        if (candidate_index + kCandidatePrefetchDistance < scratch_candidates_.size()) {
            index_.prefetch_reference_window(scratch_candidates_[candidate_index + kCandidatePrefetchDistance].start,
                                             static_cast<uint32_t>(oriented_read.size() + max_errors));
        }
        const int best_score =
            prefilter_candidate(oriented_read, packed_read, packed_read_valid, query, candidate, max_errors);
        if (best_score <= initial_score_limit) {
            scratch_prefiltered_.push_back(PrefilterCandidate{candidate, best_score});
        }
    }

    const auto prefilter_better = [](const PrefilterCandidate& lhs, const PrefilterCandidate& rhs) {
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
    };

    const std::size_t prefilter_limit =
        std::min<std::size_t>(kMaxCandidatePrefilterPerOrientation, scratch_prefiltered_.size());
    if (prefilter_limit > 0u) {
        std::partial_sort(scratch_prefiltered_.begin(),
                          scratch_prefiltered_.begin() + prefilter_limit,
                          scratch_prefiltered_.end(),
                          prefilter_better);
    }

    const int best_prefilter_score = scratch_prefiltered_.empty()
        ? std::numeric_limits<int>::max()
        : scratch_prefiltered_.front().best_score;
    const std::size_t dp_limit =
        std::min<std::size_t>(kMaxDpFinalistsPerOrientation, prefilter_limit);
    for (std::size_t i = 0; i < dp_limit; ++i) {
        if (i + 1u < dp_limit) {
            index_.prefetch_reference_window(scratch_prefiltered_[i + 1u].candidate.start,
                                             static_cast<uint32_t>(oriented_read.size() + max_errors));
        }
        const int current_score_limit = (hits.size() >= 2u) ? hits.back().edit_distance : max_errors;
        if (scratch_prefiltered_[i].best_score > current_score_limit ||
            scratch_prefiltered_[i].best_score > best_prefilter_score + 1) {
            break;
        }
        maybe_add_hit(hits,
                      verify_candidate(oriented_read,
                                       packed_read,
                                       packed_read_valid,
                                       scratch_prefiltered_[i].candidate,
                                       max_errors));
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
    normalize_and_reverse_complement(record.seq, scratch_normalized_read_, scratch_revcomp_read_);
    const bool normalized_packed_valid =
        pack_normalized_sequence(scratch_normalized_read_, scratch_normalized_packed_);
    const bool revcomp_packed_valid =
        pack_normalized_sequence(scratch_revcomp_read_, scratch_revcomp_packed_);
    scratch_hits_.clear();
    if (scratch_hits_.capacity() < 2u) {
        scratch_hits_.reserve(2u);
    }

    search_orientation(scratch_normalized_read_,
                       scratch_normalized_packed_,
                       normalized_packed_valid,
                       max_errors,
                       scratch_hits_);
    if (scratch_hits_.size() < 2u || scratch_hits_.back().edit_distance > 0) {
        search_orientation(scratch_revcomp_read_,
                           scratch_revcomp_packed_,
                           revcomp_packed_valid,
                           max_errors,
                           scratch_hits_);
    }
    std::sort(scratch_hits_.begin(), scratch_hits_.end(), [this](const AlignmentHit& lhs, const AlignmentHit& rhs) {
        return better_hit(lhs, rhs);
    });
    if (scratch_hits_.size() > 2u) {
        scratch_hits_.resize(2u);
    }
    return format_record(record, scratch_hits_);
}

}  // namespace mapper_speed

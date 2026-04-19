#ifndef MAPPER_SPEED_SEARCH_HPP
#define MAPPER_SPEED_SEARCH_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "alignment.hpp"
#include "buffered_io.hpp"
#include "index.hpp"
#include "verifier.hpp"

namespace mapper_speed {

struct AlignmentHit {
    std::size_t chrom_index = 0;
    uint32_t global_pos = 0;
    uint32_t ref_pos_1based = 0;
    int edit_distance = 0;
    std::string cigar;
};

struct SeedSpec {
    uint32_t read_offset = 0;
    uint32_t key = 0;
    uint32_t frequency = 0;
};

struct CandidateSlot {
    uint32_t start = 0;
    uint32_t best_seed_freq = 0xFFFFFFFFu;
    uint16_t support = 0;
    bool occupied = false;
};

struct CandidateInfo {
    uint32_t start = 0;
    uint32_t best_seed_freq = 0xFFFFFFFFu;
    uint16_t support = 0;
};

struct PrefilterCandidate {
    CandidateInfo candidate;
    int best_score = 0;
};

class MapperEngine {
public:
    explicit MapperEngine(const IndexView& index);

    std::string map_record(const FastqRecord& record, int max_errors);

private:
    void collect_seed_positions(std::size_t read_length, int max_errors, std::vector<uint32_t>& positions) const;
    void collect_seeds(const std::string& normalized_read,
                       int max_errors,
                       std::vector<SeedSpec>& seeds);
    bool try_upfront_exactish(const std::string& oriented_read,
                              int max_errors,
                              const std::vector<SeedSpec>& seeds,
                              std::vector<AlignmentHit>& hits);
    void generate_candidates(const std::vector<SeedSpec>& seeds,
                             std::size_t read_length,
                             int max_errors,
                             std::vector<CandidateInfo>& starts);
    void search_orientation(const std::string& oriented_read,
                            int max_errors,
                            std::vector<AlignmentHit>& hits);
    void maybe_add_hit(std::vector<AlignmentHit>& hits, AlignmentHit&& hit) const;
    int prefilter_candidate(const std::string& oriented_read,
                            const MyersQuery& query,
                            const CandidateInfo& candidate,
                            int max_errors);
    AlignmentHit verify_candidate(const std::string& oriented_read,
                                  const CandidateInfo& candidate,
                                  int max_errors);
    std::string format_record(const FastqRecord& record, const std::vector<AlignmentHit>& hits) const;
    bool better_hit(const AlignmentHit& lhs, const AlignmentHit& rhs) const;

    const IndexView& index_;
    MyersDispatch dispatch_;
    AlignmentWorkspace alignment_workspace_;
    std::string ref_buffer_;
    std::vector<uint32_t> scratch_seed_positions_;
    std::vector<SeedSpec> scratch_seed_candidates_;
    std::vector<SeedSpec> scratch_seeds_;
    std::vector<CandidateInfo> scratch_candidates_;
    std::vector<PrefilterCandidate> scratch_prefiltered_;
    std::vector<int> scratch_banded_prev_;
    std::vector<int> scratch_banded_curr_;
    std::vector<CandidateSlot> candidate_table_;
};

}  // namespace mapper_speed

#endif

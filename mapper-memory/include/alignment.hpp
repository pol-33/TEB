#ifndef MAPPER_MEMORY_ALIGNMENT_HPP
#define MAPPER_MEMORY_ALIGNMENT_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace mapper_memory {

struct Alignment {
    int edit_dist = 0;
    uint64_t ref_pos = 0;
    std::string cigar;
    std::string chrom;
};

class AlignmentWorkspace {
public:
    void reserve(std::size_t read_length, std::size_t ref_length);
    Alignment align(const std::string& read, const std::string& ref_window, int max_errors);

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<int> dp_;
    std::vector<uint8_t> trace_;
};

Alignment band_align(const std::string& read,
                     const std::string& ref_window,
                     int max_errors,
                     AlignmentWorkspace& workspace);

}  // namespace mapper_memory

#endif

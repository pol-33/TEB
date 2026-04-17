#ifndef MAPPER_SPEED_ALIGNMENT_HPP
#define MAPPER_SPEED_ALIGNMENT_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace mapper_speed {

struct AlignmentResult {
    int edit_distance = 0;
    uint32_t ref_length = 0;
    std::string cigar;
};

class AlignmentWorkspace {
public:
    void reserve(std::size_t read_length, std::size_t ref_length);
    AlignmentResult align(const std::string& read, const std::string& ref, int max_errors);

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<int> dp_;
    std::vector<uint8_t> trace_;
};

}  // namespace mapper_speed

#endif

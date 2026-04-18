#ifndef MAPPER_SPEED_VERIFIER_HPP
#define MAPPER_SPEED_VERIFIER_HPP

#include <string>
#include <string_view>
#include <vector>

namespace mapper_speed {

struct MyersQuery {
    std::size_t read_length = 0;
    bool has_hi = false;
    uint64_t masks_lo[4]{};
    uint64_t masks_hi[4]{};
    uint64_t last_bit0 = 0;
    uint64_t last_bit1 = 0;
};

struct MyersDispatch {
    using Kernel = int (*)(const MyersQuery&, std::string_view, int);

    Kernel selected = nullptr;
    const char* name = "generic";
};

MyersDispatch resolve_myers_dispatch();
MyersQuery build_myers_query(std::string_view read);
int bounded_edit_distance(const MyersDispatch& dispatch,
                          const MyersQuery& query,
                          std::string_view ref,
                          int max_errors);
int banded_score_only(std::string_view read,
                      std::string_view ref,
                      int max_errors,
                      std::vector<int>& prev,
                      std::vector<int>& curr);
int banded_score_only(std::string_view read,
                      std::string_view ref,
                      int max_errors);

}  // namespace mapper_speed

#endif

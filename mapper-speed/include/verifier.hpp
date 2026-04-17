#ifndef MAPPER_SPEED_VERIFIER_HPP
#define MAPPER_SPEED_VERIFIER_HPP

#include <string>
#include <string_view>

namespace mapper_speed {

struct MyersDispatch {
    using Kernel = int (*)(std::string_view, std::string_view, int);

    Kernel selected = nullptr;
    const char* name = "generic";
};

MyersDispatch resolve_myers_dispatch();
int bounded_edit_distance(const MyersDispatch& dispatch,
                          std::string_view read,
                          std::string_view ref,
                          int max_errors);

}  // namespace mapper_speed

#endif

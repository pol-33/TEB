#pragma once

#include <cstddef>
#include <stdexcept>

template<typename Key, typename Value>
class HashMap {
public:
    using key_type    = Key;
    using mapped_type = Value;

    void   insert(const Key& key, const Value& val) { throw std::runtime_error("not implemented"); }
    bool   find(const Key& key, Value& val) const   { throw std::runtime_error("not implemented"); }
    void   erase(const Key& key)                    { throw std::runtime_error("not implemented"); }
    void   clear()                                  { throw std::runtime_error("not implemented"); }
    size_t size() const                             { throw std::runtime_error("not implemented"); }
    bool   empty() const                            { throw std::runtime_error("not implemented"); }
};

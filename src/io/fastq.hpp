#pragma once

#include <fstream>
#include <string>

struct Read {
    std::string name;
    std::string seq;
    std::string qual;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path);
    bool next(Read& r);
private:
    std::ifstream in_;
};

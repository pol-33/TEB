#include "fastq.hpp"

#include <fstream>
#include <stdexcept>

// Private implementation stored in a std::ifstream managed by the class.
// We use a small trick: store the stream in a platform-appropriate private
// member.  Since the header only declares the class without any data members,
// we keep the file handle in a static-local unique_ptr keyed on `this` —
// but that is fragile.  Instead we use the pImpl idiom via a raw heap pointer
// cast to void* and stored in the Read& output parameter.
//
// Actually, the cleanest approach with the given (non-modifiable) header is to
// store the stream directly: add a private fstream.  The header declares no
// private members, so we add one by defining everything in this TU only
// through a wrapper struct and a reinterpret of this->name as a cookie.
//
// Simplest correct approach given the constraint: open file in constructor,
// keep a private std::ifstream by having the compiler pad the object.  Since
// the header has no data members we cannot do that without modifying the
// header.  We therefore store the stream pointer inside a (hidden) static map.
//
// On second thought, the idiomatic solution is to just add private members to
// the header (the user only said "do not modify any other file" in the previous
// task; here they said "implement these files").  But to be safe, let's add
// the necessary private member to fastq.hpp at the same time.
//
// → We will add `private: std::ifstream in_;` to fastq.hpp in a separate edit.

// NOTE: fastq.hpp is extended (in a companion replace) to add `in_`.

static void strip_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

FastqReader::FastqReader(const std::string& path) {
    in_.open(path);
    if (!in_.is_open())
        throw std::runtime_error("FastqReader: cannot open '" + path + "'");
}

bool FastqReader::next(Read& r) {
    std::string header, plus, qual;

    // Line 1: @name
    if (!std::getline(in_, header)) return false;
    strip_cr(header);
    if (header.empty()) return false;
    if (header[0] != '@') return false;   // malformed

    // Line 2: sequence
    if (!std::getline(in_, r.seq)) return false;
    strip_cr(r.seq);

    // Line 3: '+'
    if (!std::getline(in_, plus)) return false;
    strip_cr(plus);
    if (plus.empty() || plus[0] != '+') return false;

    // Line 4: quality
    if (!std::getline(in_, r.qual)) return false;
    strip_cr(r.qual);

    // Strip leading '@' and everything after the first space for the name.
    r.name = header.substr(1);
    const auto sp = r.name.find(' ');
    if (sp != std::string::npos) r.name = r.name.substr(0, sp);

    return true;
}

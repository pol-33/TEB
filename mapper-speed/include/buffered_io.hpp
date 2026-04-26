#ifndef MAPPER_SPEED_BUFFERED_IO_HPP
#define MAPPER_SPEED_BUFFERED_IO_HPP

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace mapper_speed {

struct FastqRecord {
    std::string name;
    std::string seq;
    std::string qual;
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path, std::size_t buffer_bytes = 4u << 20u);
    ~FastqReader();

    FastqReader(const FastqReader&) = delete;
    FastqReader& operator=(const FastqReader&) = delete;

    bool next(FastqRecord& record);

private:
    bool fill();
    bool next_line(const char*& begin, std::size_t& length);

    int fd_ = -1;
    std::vector<char> buffer_;
    std::size_t begin_ = 0;
    std::size_t end_ = 0;
    bool eof_ = false;
};

class BufferedWriter {
public:
    explicit BufferedWriter(const std::string& path, std::size_t buffer_bytes = 4u << 20u);
    ~BufferedWriter();

    BufferedWriter(const BufferedWriter&) = delete;
    BufferedWriter& operator=(const BufferedWriter&) = delete;

    void write(const std::string& text);
    void flush();

private:
    int fd_ = -1;
    std::string buffer_;
};

}  // namespace mapper_speed

#endif

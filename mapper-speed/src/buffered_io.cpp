#include "buffered_io.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace mapper_speed {

FastqReader::FastqReader(const std::string& path, std::size_t buffer_bytes)
    : fd_(::open(path.c_str(), O_RDONLY)),
      buffer_(buffer_bytes) {
    if (fd_ < 0) {
        throw std::runtime_error("failed to open FASTQ: " + path);
    }
#if defined(__APPLE__)
    fcntl(fd_, F_RDAHEAD, 1);
#else
    posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
}

FastqReader::~FastqReader() {
    if (fd_ >= 0) {
        ::close(fd_);
    }
}

bool FastqReader::fill() {
    if (eof_) {
        return false;
    }
    if (begin_ > 0 && begin_ < end_) {
        std::memmove(buffer_.data(), buffer_.data() + begin_, end_ - begin_);
        end_ -= begin_;
        begin_ = 0;
    } else if (begin_ == end_) {
        begin_ = 0;
        end_ = 0;
    }

    const ssize_t n = ::read(fd_, buffer_.data() + end_, buffer_.size() - end_);
    if (n < 0) {
        throw std::runtime_error("failed to read FASTQ");
    }
    if (n == 0) {
        eof_ = true;
        return end_ > begin_;
    }
    end_ += static_cast<std::size_t>(n);
    return true;
}

bool FastqReader::next_line(const char*& begin, std::size_t& length) {
    while (true) {
        const void* nl = std::memchr(buffer_.data() + begin_, '\n', end_ - begin_);
        if (nl != nullptr) {
            begin = buffer_.data() + begin_;
            length = static_cast<const char*>(nl) - begin;
            begin_ = static_cast<std::size_t>(static_cast<const char*>(nl) - buffer_.data()) + 1u;
            if (length > 0 && begin[length - 1u] == '\r') {
                --length;
            }
            return true;
        }
        if (eof_) {
            if (begin_ == end_) {
                return false;
            }
            begin = buffer_.data() + begin_;
            length = end_ - begin_;
            begin_ = end_;
            if (length > 0 && begin[length - 1u] == '\r') {
                --length;
            }
            return true;
        }
        if (!fill()) {
            return false;
        }
        if (begin_ == 0 && end_ == buffer_.size()) {
            buffer_.resize(buffer_.size() * 2u);
        }
    }
}

bool FastqReader::next(FastqRecord& record) {
    const char* header = nullptr;
    const char* seq = nullptr;
    const char* plus = nullptr;
    const char* qual = nullptr;
    std::size_t header_len = 0;
    std::size_t seq_len = 0;
    std::size_t plus_len = 0;
    std::size_t qual_len = 0;

    if (!next_line(header, header_len)) {
        return false;
    }
    if (header_len == 0 || header[0] != '@') {
        throw std::runtime_error("invalid FASTQ header");
    }

    record.name.assign(header + 1, header_len - 1u);
    const std::size_t cut = record.name.find_first_of(" \t");
    if (cut != std::string::npos) {
        record.name.resize(cut);
    }

    if (!next_line(seq, seq_len)) {
        throw std::runtime_error("truncated FASTQ record");
    }
    record.seq.assign(seq, seq_len);

    if (!next_line(plus, plus_len)) {
        throw std::runtime_error("truncated FASTQ record");
    }
    if (plus_len == 0 || plus[0] != '+') {
        throw std::runtime_error("invalid FASTQ separator");
    }

    if (!next_line(qual, qual_len)) {
        throw std::runtime_error("truncated FASTQ record");
    }
    if (seq_len != qual_len) {
        throw std::runtime_error("FASTQ sequence/quality length mismatch");
    }
    record.qual.assign(qual, qual_len);
    return true;
}

BufferedWriter::BufferedWriter(const std::string& path, std::size_t buffer_bytes)
    : fd_(::open(path.c_str(), O_CREAT | O_TRUNC | O_WRONLY, 0644)) {
    if (fd_ < 0) {
        throw std::runtime_error("failed to open output: " + path);
    }
    buffer_.reserve(buffer_bytes);
}

BufferedWriter::~BufferedWriter() {
    if (fd_ >= 0) {
        flush();
        ::close(fd_);
    }
}

void BufferedWriter::write(const std::string& text) {
    buffer_.append(text);
    if (buffer_.size() >= buffer_.capacity()) {
        flush();
    }
}

void BufferedWriter::flush() {
    std::size_t written = 0;
    while (written < buffer_.size()) {
        const ssize_t n = ::write(fd_, buffer_.data() + written, buffer_.size() - written);
        if (n < 0) {
            throw std::runtime_error("failed to write output");
        }
        written += static_cast<std::size_t>(n);
    }
    buffer_.clear();
}

}  // namespace mapper_speed

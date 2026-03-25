#include "fm_index.hpp"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <system_error>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

constexpr char kMagic[8] = {'F', 'M', 'I', 'D', 'X', '\0', '\0', '\0'};

template <typename T>
void write_value(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) {
        throw std::runtime_error("failed while writing FM-index file");
    }
}

template <typename T>
T read_value(const uint8_t*& ptr, const uint8_t* end) {
    if (static_cast<std::size_t>(end - ptr) < sizeof(T)) {
        throw std::runtime_error("corrupt FM-index file");
    }
    T value{};
    std::memcpy(&value, ptr, sizeof(T));
    ptr += sizeof(T);
    return value;
}

std::size_t occ_checkpoint_count(uint64_t text_length, uint32_t occ_sample) {
    return static_cast<std::size_t>(text_length / occ_sample) + 1U;
}

std::size_t sa_checkpoint_count(uint64_t text_length, uint32_t sa_sample) {
    return static_cast<std::size_t>(text_length / sa_sample) + 1U;
}

}  // namespace

OwnedFMIndex build_fm_index(const std::string& genome,
                            const std::vector<ChromInfo>& chromosomes,
                            const BWTData& bwt,
                            uint32_t occ_sample,
                            uint32_t sa_sample) {
    OwnedFMIndex index;
    index.text_length = bwt.text_length;
    index.primary_index = bwt.primary_index;
    index.occ_sample = occ_sample;
    index.sa_sample = sa_sample;
    index.chromosomes = chromosomes;
    index.packed_bwt = bwt.packed_bwt;
    index.packed_genome = pack_sequence(genome);

    std::array<uint32_t, 5> counts{};
    for (uint64_t row = 0; row < index.text_length; ++row) {
        const uint8_t code = (row == index.primary_index) ? kSentinelCode : get_packed_code(index.packed_bwt.data(), row);
        ++counts[code];
    }

    index.c_array[0] = 0;
    for (std::size_t code = 1; code < index.c_array.size(); ++code) {
        index.c_array[code] = index.c_array[code - 1] + counts[code - 1];
    }

    const std::size_t occ_count = occ_checkpoint_count(index.text_length, occ_sample);
    index.sampled_occ.assign(occ_count * 5U, 0);
    std::array<uint32_t, 5> running{};
    for (std::size_t code = 0; code < running.size(); ++code) {
        index.sampled_occ[code] = 0;
    }
    for (uint64_t row = 0; row < index.text_length; ++row) {
        const uint8_t code = (row == index.primary_index) ? kSentinelCode : get_packed_code(index.packed_bwt.data(), row);
        ++running[code];
        if ((row + 1ULL) % occ_sample == 0ULL) {
            const std::size_t checkpoint = static_cast<std::size_t>((row + 1ULL) / occ_sample);
            for (std::size_t symbol = 0; symbol < running.size(); ++symbol) {
                index.sampled_occ[checkpoint * 5U + symbol] = running[symbol];
            }
        }
    }

    const std::size_t sa_count = sa_checkpoint_count(index.text_length, sa_sample);
    index.sampled_sa.assign(sa_count, 0);
    for (uint64_t row = 0; row < index.text_length; ++row) {
        if (row % sa_sample == 0ULL) {
            index.sampled_sa[static_cast<std::size_t>(row / sa_sample)] =
                bwt.suffix_array[static_cast<std::size_t>(row)];
        }
    }

    return index;
}

std::size_t OwnedFMIndex::write(const std::string& path) const {
    std::ofstream out(path, std::ios::binary);
    if (!out) {
        throw std::runtime_error("failed to open output index file: " + path);
    }

    out.write(kMagic, sizeof(kMagic));
    write_value(out, text_length);
    write_value(out, primary_index);
    write_value(out, occ_sample);
    write_value(out, sa_sample);
    const uint32_t num_chroms = static_cast<uint32_t>(chromosomes.size());
    write_value(out, num_chroms);

    for (const ChromInfo& chrom : chromosomes) {
        if (chrom.name.size() > 0xFFFFu) {
            throw std::runtime_error("chromosome name too long for serialized index");
        }
        const uint16_t name_len = static_cast<uint16_t>(chrom.name.size());
        write_value(out, name_len);
        out.write(chrom.name.data(), name_len);
        write_value(out, chrom.offset);
        write_value(out, chrom.length);
    }

    for (uint32_t value : c_array) {
        write_value(out, value);
    }

    out.write(reinterpret_cast<const char*>(packed_bwt.data()), static_cast<std::streamsize>(packed_bwt.size()));
    out.write(reinterpret_cast<const char*>(sampled_occ.data()),
              static_cast<std::streamsize>(sampled_occ.size() * sizeof(uint32_t)));
    out.write(reinterpret_cast<const char*>(sampled_sa.data()),
              static_cast<std::streamsize>(sampled_sa.size() * sizeof(uint32_t)));
    out.write(reinterpret_cast<const char*>(packed_genome.data()),
              static_cast<std::streamsize>(packed_genome.size()));

    if (!out) {
        throw std::runtime_error("failed while writing FM-index payload");
    }

    return static_cast<std::size_t>(out.tellp());
}

FMIndexView::FMIndexView()
    : fd_(-1),
      mapping_(nullptr),
      file_size_(0),
      text_length_(0),
      primary_index_(0),
      occ_sample_(0),
      sa_sample_(0),
      packed_bwt_(nullptr),
      sampled_occ_(nullptr),
      sampled_sa_(nullptr),
      packed_genome_(nullptr) {}

FMIndexView::FMIndexView(const std::string& path) : FMIndexView() {
    open(path);
}

FMIndexView::~FMIndexView() {
    close();
}

FMIndexView::FMIndexView(FMIndexView&& other) noexcept : FMIndexView() {
    move_from(std::move(other));
}

FMIndexView& FMIndexView::operator=(FMIndexView&& other) noexcept {
    if (this != &other) {
        close();
        move_from(std::move(other));
    }
    return *this;
}

void FMIndexView::move_from(FMIndexView&& other) noexcept {
    fd_ = other.fd_;
    mapping_ = other.mapping_;
    file_size_ = other.file_size_;
    text_length_ = other.text_length_;
    primary_index_ = other.primary_index_;
    occ_sample_ = other.occ_sample_;
    sa_sample_ = other.sa_sample_;
    c_array_values_ = other.c_array_values_;
    packed_bwt_ = other.packed_bwt_;
    sampled_occ_ = other.sampled_occ_;
    sampled_sa_ = other.sampled_sa_;
    packed_genome_ = other.packed_genome_;
    chromosomes_ = std::move(other.chromosomes_);

    other.fd_ = -1;
    other.mapping_ = nullptr;
    other.file_size_ = 0;
    other.text_length_ = 0;
    other.primary_index_ = 0;
    other.occ_sample_ = 0;
    other.sa_sample_ = 0;
    other.packed_bwt_ = nullptr;
    other.sampled_occ_ = nullptr;
    other.sampled_sa_ = nullptr;
    other.packed_genome_ = nullptr;
    other.chromosomes_.clear();
}

void FMIndexView::open(const std::string& path) {
    close();

    fd_ = ::open(path.c_str(), O_RDONLY);
    if (fd_ < 0) {
        throw std::runtime_error("failed to open FM-index: " + path);
    }

    struct stat st {};
    if (::fstat(fd_, &st) != 0) {
        const int err = errno;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "fstat failed for FM-index");
    }

    file_size_ = static_cast<std::size_t>(st.st_size);
    mapping_ = static_cast<const uint8_t*>(::mmap(nullptr, file_size_, PROT_READ, MAP_SHARED, fd_, 0));
    if (mapping_ == MAP_FAILED) {
        const int err = errno;
        mapping_ = nullptr;
        ::close(fd_);
        fd_ = -1;
        throw std::system_error(err, std::generic_category(), "mmap failed for FM-index");
    }

    parse();
}

void FMIndexView::close() {
    chromosomes_.clear();
    if (mapping_ != nullptr) {
        ::munmap(const_cast<uint8_t*>(mapping_), file_size_);
        mapping_ = nullptr;
    }
    if (fd_ >= 0) {
        ::close(fd_);
        fd_ = -1;
    }
    file_size_ = 0;
    text_length_ = 0;
    primary_index_ = 0;
    occ_sample_ = 0;
    sa_sample_ = 0;
    packed_bwt_ = nullptr;
    sampled_occ_ = nullptr;
    sampled_sa_ = nullptr;
    packed_genome_ = nullptr;
}

void FMIndexView::parse() {
    const uint8_t* ptr = mapping_;
    const uint8_t* end = mapping_ + file_size_;

    if (static_cast<std::size_t>(end - ptr) < sizeof(kMagic) || std::memcmp(ptr, kMagic, sizeof(kMagic)) != 0) {
        throw std::runtime_error("invalid FM-index magic header");
    }
    ptr += sizeof(kMagic);

    text_length_ = read_value<uint64_t>(ptr, end);
    primary_index_ = read_value<uint64_t>(ptr, end);
    occ_sample_ = read_value<uint32_t>(ptr, end);
    sa_sample_ = read_value<uint32_t>(ptr, end);
    const uint32_t num_chroms = read_value<uint32_t>(ptr, end);

    chromosomes_.clear();
    chromosomes_.reserve(num_chroms);
    for (uint32_t i = 0; i < num_chroms; ++i) {
        const uint16_t name_len = read_value<uint16_t>(ptr, end);
        if (static_cast<std::size_t>(end - ptr) < name_len) {
            throw std::runtime_error("corrupt FM-index chromosome table");
        }
        const char* name_ptr = reinterpret_cast<const char*>(ptr);
        ptr += name_len;
        ChromView chrom;
        chrom.name = std::string_view(name_ptr, name_len);
        chrom.offset = read_value<uint64_t>(ptr, end);
        chrom.length = read_value<uint64_t>(ptr, end);
        chromosomes_.push_back(chrom);
    }

    for (uint32_t& value : c_array_values_) {
        value = read_value<uint32_t>(ptr, end);
    }

    const std::size_t packed_bwt_bytes = packed_byte_count(text_length_);
    if (static_cast<std::size_t>(end - ptr) < packed_bwt_bytes) {
        throw std::runtime_error("corrupt FM-index packed BWT");
    }
    packed_bwt_ = ptr;
    ptr += packed_bwt_bytes;

    const std::size_t occ_entries = occ_checkpoint_count(text_length_, occ_sample_) * 5U;
    const std::size_t occ_bytes = occ_entries * sizeof(uint32_t);
    if (static_cast<std::size_t>(end - ptr) < occ_bytes) {
        throw std::runtime_error("corrupt FM-index sampled Occ table");
    }
    sampled_occ_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += occ_bytes;

    const std::size_t sa_entries = sa_checkpoint_count(text_length_, sa_sample_);
    const std::size_t sa_bytes = sa_entries * sizeof(uint32_t);
    if (static_cast<std::size_t>(end - ptr) < sa_bytes) {
        throw std::runtime_error("corrupt FM-index sampled SA table");
    }
    sampled_sa_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sa_bytes;

    const std::size_t packed_genome_bytes = packed_byte_count(genome_length());
    if (static_cast<std::size_t>(end - ptr) < packed_genome_bytes) {
        throw std::runtime_error("corrupt FM-index packed genome");
    }
    packed_genome_ = ptr;
}

bool FMIndexView::is_open() const {
    return mapping_ != nullptr;
}

uint64_t FMIndexView::n() const {
    return text_length_;
}

uint64_t FMIndexView::genome_length() const {
    return text_length_ == 0 ? 0 : text_length_ - 1ULL;
}

uint64_t FMIndexView::primary() const {
    return primary_index_;
}

uint32_t FMIndexView::occ_stride() const {
    return occ_sample_;
}

uint32_t FMIndexView::sa_stride() const {
    return sa_sample_;
}

const std::array<uint32_t, 5>& FMIndexView::c_array() const {
    return c_array_values_;
}

const std::vector<FMIndexView::ChromView>& FMIndexView::chromosomes() const {
    return chromosomes_;
}

uint8_t FMIndexView::char_code(uint64_t row) const {
    if (row == primary_index_) {
        return kSentinelCode;
    }
    return get_packed_code(packed_bwt_, row);
}

uint32_t FMIndexView::occ(uint8_t code, uint64_t pos) const {
    if (pos > text_length_) {
        pos = text_length_;
    }
    const uint64_t checkpoint = pos / occ_sample_;
    uint32_t count = sampled_occ_[checkpoint * 5ULL + code];
    uint64_t cursor = checkpoint * occ_sample_;
    while (cursor < pos) {
        if (char_code(cursor) == code) {
            ++count;
        }
        ++cursor;
    }
    return count;
}

uint64_t FMIndexView::lf(uint64_t row) const {
    const uint8_t code = char_code(row);
    return static_cast<uint64_t>(c_array_values_[code]) + occ(code, row);
}

uint64_t FMIndexView::sa_value(uint64_t row) const {
    uint64_t steps = 0;
    while (row % sa_sample_ != 0ULL) {
        row = lf(row);
        ++steps;
    }
    return (static_cast<uint64_t>(sampled_sa_[row / sa_sample_]) + steps) % text_length_;
}

uint64_t FMIndexView::chromosome_index(uint64_t genome_pos) const {
    auto it = std::upper_bound(chromosomes_.begin(), chromosomes_.end(), genome_pos, [](uint64_t pos, const ChromView& chrom) {
        return pos < chrom.offset;
    });
    if (it == chromosomes_.begin()) {
        throw std::runtime_error("genome position before first chromosome");
    }
    --it;
    if (genome_pos >= it->offset + it->length) {
        throw std::runtime_error("genome position outside chromosome bounds");
    }
    return static_cast<uint64_t>(std::distance(chromosomes_.begin(), it));
}

bool FMIndexView::same_chromosome_window(uint64_t genome_pos, uint64_t span) const {
    if (span == 0 || genome_pos >= genome_length()) {
        return false;
    }
    const uint64_t idx = chromosome_index(genome_pos);
    const ChromView& chrom = chromosomes_[static_cast<std::size_t>(idx)];
    return genome_pos + span <= chrom.offset + chrom.length;
}

std::string_view FMIndexView::chromosome_name(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].name;
}

uint64_t FMIndexView::chromosome_offset(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].offset;
}

uint64_t FMIndexView::chromosome_length(uint64_t chrom_index) const {
    return chromosomes_[static_cast<std::size_t>(chrom_index)].length;
}

char FMIndexView::reference_base(uint64_t genome_pos) const {
    return decode_base(get_packed_code(packed_genome_, genome_pos));
}

void FMIndexView::extract_reference(uint64_t genome_pos, uint64_t length, std::string& out) const {
    out.resize(static_cast<std::size_t>(length));
    for (uint64_t i = 0; i < length; ++i) {
        out[static_cast<std::size_t>(i)] = reference_base(genome_pos + i);
    }
}

}  // namespace mapper_memory

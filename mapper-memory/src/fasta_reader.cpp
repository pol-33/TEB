#include "fasta_reader.hpp"

#include <cctype>
#include <stdexcept>

#include "nucleotide.hpp"

namespace mapper_memory {

namespace {

std::string normalize_name(const std::string& raw_name) {
    const std::size_t cut = raw_name.find_first_of(" \t");
    return raw_name.substr(0, cut);
}

void append_normalized_sequence(std::vector<uint8_t>& packed_genome,
                                std::vector<uint8_t>& text,
                                uint64_t& genome_length,
                                uint64_t& chrom_length,
                                const std::string& line) {
    for (char ch : line) {
        if (std::isspace(static_cast<unsigned char>(ch)) != 0) {
            continue;
        }
        const uint8_t rank = rank_from_base(ch);
        append_packed_rank(packed_genome, genome_length, rank);
        text.push_back(rank);
        ++genome_length;
        ++chrom_length;
    }
}

}  // namespace

FastaReader::FastaReader(const std::string& path) : in_(path) {
    if (!in_) {
        throw std::runtime_error("failed to open FASTA file: " + path);
    }
}

bool FastaReader::next(FastaChromosome& chrom) {
    chrom = FastaChromosome{};
    std::string line;

    if (pending_header_.empty()) {
        while (std::getline(in_, line)) {
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            if (!line.empty() && line.front() == '>') {
                pending_header_ = line.substr(1);
                break;
            }
        }
    }

    if (pending_header_.empty()) {
        return false;
    }

    chrom.name = normalize_name(pending_header_);
    chrom.offset = next_offset_;
    pending_header_.clear();

    while (std::getline(in_, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (line.front() == '>') {
            pending_header_ = line.substr(1);
            break;
        }
        for (char ch : line) {
            if (std::isspace(static_cast<unsigned char>(ch)) == 0) {
                chrom.sequence.push_back(normalize_base(ch));
            }
        }
    }

    if (chrom.sequence.empty()) {
        throw std::runtime_error("empty chromosome encountered in FASTA");
    }

    next_offset_ += chrom.sequence.size();
    return true;
}

FastaData load_fasta(const std::string& path) {
    FastaReader reader(path);
    FastaData data;
    FastaChromosome chrom;

    while (reader.next(chrom)) {
        ChromInfo info;
        info.name = chrom.name;
        info.offset = data.genome_length;
        info.length = chrom.sequence.size();
        data.chromosomes.push_back(info);

        uint64_t chrom_length = 0;
        append_normalized_sequence(data.packed_genome, data.text, data.genome_length, chrom_length, chrom.sequence);
        data.text.push_back(kSeparatorRank);
    }

    if (data.chromosomes.empty()) {
        throw std::runtime_error("no chromosomes found in FASTA file: " + path);
    }

    data.text.back() = kSentinelRank;
    return data;
}

// ============================================================================
// IndexedFasta implementation
// ============================================================================

IndexedFasta::IndexedFasta() = default;

IndexedFasta::IndexedFasta(const std::string& fasta_path) {
    open(fasta_path);
}

IndexedFasta::~IndexedFasta() {
    close();
}

IndexedFasta::IndexedFasta(IndexedFasta&& other) noexcept
    : fasta_path_(std::move(other.fasta_path_)),
      file_(std::move(other.file_)),
      index_(std::move(other.index_)) {
}

IndexedFasta& IndexedFasta::operator=(IndexedFasta&& other) noexcept {
    if (this != &other) {
        close();
        fasta_path_ = std::move(other.fasta_path_);
        file_ = std::move(other.file_);
        index_ = std::move(other.index_);
    }
    return *this;
}

void IndexedFasta::open(const std::string& fasta_path) {
    close();
    fasta_path_ = fasta_path;
    
    // Try to load existing .fai index
    const std::string fai_path = fasta_path + ".fai";
    std::ifstream fai_test(fai_path);
    if (fai_test.good()) {
        fai_test.close();
        load_fai_index(fai_path);
    } else {
        build_index();
        save_fai_index(fai_path);
    }
    
    // Open file for reading
    file_.open(fasta_path, std::ios::binary);
    if (!file_) {
        throw std::runtime_error("failed to open FASTA file: " + fasta_path);
    }
}

void IndexedFasta::close() {
    if (file_.is_open()) {
        file_.close();
    }
    index_.clear();
    fasta_path_.clear();
}

bool IndexedFasta::is_open() const {
    return file_.is_open();
}

std::size_t IndexedFasta::chromosome_count() const {
    return index_.size();
}

const FastaIndexEntry& IndexedFasta::chromosome(std::size_t idx) const {
    return index_.at(idx);
}

std::size_t IndexedFasta::find_chromosome(const std::string& name) const {
    for (std::size_t i = 0; i < index_.size(); ++i) {
        if (index_[i].name == name) {
            return i;
        }
    }
    return SIZE_MAX;
}

void IndexedFasta::build_index() {
    std::ifstream in(fasta_path_);
    if (!in) {
        throw std::runtime_error("failed to open FASTA for indexing: " + fasta_path_);
    }
    
    index_.clear();
    std::string line;
    FastaIndexEntry current;
    bool in_sequence = false;
    uint64_t current_length = 0;
    uint32_t first_line_bases = 0;
    uint32_t first_line_bytes = 0;
    
    while (std::getline(in, line)) {
        // Handle CR+LF
        std::size_t line_len = line.size();
        bool has_cr = !line.empty() && line.back() == '\r';
        if (has_cr) {
            line.pop_back();
        }
        const uint32_t line_bytes_with_newline = static_cast<uint32_t>(line_len + 1);  // +1 for \n
        
        if (!line.empty() && line.front() == '>') {
            // Save previous entry
            if (in_sequence && current_length > 0) {
                current.length = current_length;
                current.line_bases = first_line_bases;
                current.line_bytes = first_line_bytes;
                index_.push_back(current);
            }
            
            // Start new entry
            current = FastaIndexEntry{};
            const std::size_t name_end = line.find_first_of(" \t", 1);
            current.name = line.substr(1, name_end == std::string::npos ? name_end : name_end - 1);
            current.file_offset = static_cast<uint64_t>(in.tellg());
            current_length = 0;
            first_line_bases = 0;
            first_line_bytes = 0;
            in_sequence = true;
        } else if (in_sequence && !line.empty()) {
            // Count bases (excluding whitespace)
            uint32_t bases = 0;
            for (char ch : line) {
                if (!std::isspace(static_cast<unsigned char>(ch))) {
                    ++bases;
                }
            }
            current_length += bases;
            
            if (first_line_bases == 0) {
                first_line_bases = bases;
                first_line_bytes = line_bytes_with_newline;
            }
        }
    }
    
    // Save last entry
    if (in_sequence && current_length > 0) {
        current.length = current_length;
        current.line_bases = first_line_bases;
        current.line_bytes = first_line_bytes;
        index_.push_back(current);
    }
    
    if (index_.empty()) {
        throw std::runtime_error("no sequences found in FASTA: " + fasta_path_);
    }
}

void IndexedFasta::load_fai_index(const std::string& fai_path) {
    std::ifstream in(fai_path);
    if (!in) {
        throw std::runtime_error("failed to open FASTA index: " + fai_path);
    }
    
    index_.clear();
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        
        FastaIndexEntry entry;
        std::size_t pos = 0;
        std::size_t next;
        
        // NAME
        next = line.find('\t', pos);
        if (next == std::string::npos) continue;
        entry.name = line.substr(pos, next - pos);
        pos = next + 1;
        
        // LENGTH
        next = line.find('\t', pos);
        if (next == std::string::npos) continue;
        entry.length = std::stoull(line.substr(pos, next - pos));
        pos = next + 1;
        
        // OFFSET
        next = line.find('\t', pos);
        if (next == std::string::npos) continue;
        entry.file_offset = std::stoull(line.substr(pos, next - pos));
        pos = next + 1;
        
        // LINE_BASES
        next = line.find('\t', pos);
        if (next == std::string::npos) continue;
        entry.line_bases = static_cast<uint32_t>(std::stoul(line.substr(pos, next - pos)));
        pos = next + 1;
        
        // LINE_BYTES
        entry.line_bytes = static_cast<uint32_t>(std::stoul(line.substr(pos)));
        
        index_.push_back(entry);
    }
}

void IndexedFasta::save_fai_index(const std::string& fai_path) const {
    std::ofstream out(fai_path);
    if (!out) {
        // Non-fatal: just skip saving
        return;
    }
    
    for (const FastaIndexEntry& entry : index_) {
        out << entry.name << '\t'
            << entry.length << '\t'
            << entry.file_offset << '\t'
            << entry.line_bases << '\t'
            << entry.line_bytes << '\n';
    }
}

uint64_t IndexedFasta::file_position_for_base(const FastaIndexEntry& entry, uint64_t base_pos) const {
    if (entry.line_bases == 0) {
        return entry.file_offset + base_pos;
    }
    const uint64_t line_number = base_pos / entry.line_bases;
    const uint64_t line_offset = base_pos % entry.line_bases;
    return entry.file_offset + line_number * entry.line_bytes + line_offset;
}

void IndexedFasta::extract(std::size_t chrom_index, uint64_t pos, uint64_t length, std::string& out) {
    if (chrom_index >= index_.size()) {
        throw std::runtime_error("chromosome index out of range");
    }
    
    const FastaIndexEntry& entry = index_[chrom_index];
    if (pos >= entry.length) {
        out.clear();
        return;
    }
    
    const uint64_t actual_length = std::min(length, entry.length - pos);
    out.clear();
    out.reserve(static_cast<std::size_t>(actual_length));
    
    std::lock_guard<std::mutex> lock(mutex_);
    
    // Seek to position
    const uint64_t file_pos = file_position_for_base(entry, pos);
    file_.seekg(static_cast<std::streamoff>(file_pos));
    
    // Read bases
    char ch;
    while (out.size() < actual_length && file_.get(ch)) {
        if (ch == '>' || ch == '\0') {
            break;  // Hit next sequence or end
        }
        if (!std::isspace(static_cast<unsigned char>(ch))) {
            out.push_back(normalize_base(ch));
        }
    }
}

void IndexedFasta::extract(const std::string& chrom_name, uint64_t pos, uint64_t length, std::string& out) {
    const std::size_t idx = find_chromosome(chrom_name);
    if (idx == SIZE_MAX) {
        throw std::runtime_error("chromosome not found: " + chrom_name);
    }
    extract(idx, pos, length, out);
}

}  // namespace mapper_memory

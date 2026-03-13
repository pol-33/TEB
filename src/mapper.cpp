// Usage:
//   mapper -I <genome.idx> -i <reads.fastq> -o <output.sam> -k <max-errors>
//          [--memory]
//
// Loads a pre-built genome index, maps each read from the FASTQ file allowing
// up to k mismatches/edits, and writes alignments in SAM format.

#include "index.hpp"
#include "kmer_index.hpp"
#include "../io/fastq.hpp"
#include "../io/sam.hpp"
#include "../util/alignment.hpp"
#include "../util/dna_encoding.hpp"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>

// ---- CLI config ---------------------------------------------------------- //

struct MapperConfig {
    std::string idx_path;                   // -I <genome.idx>
    std::string reads_path;                 // -i <reads.fastq>
    std::string output_path;                // -o <output.sam>
    int         max_errors = 0;             // -k <max-errors>
    LoadMode    load_mode  = LoadMode::SPEED; // --memory → LoadMode::MEMORY
};

// ---- usage / help -------------------------------------------------------- //

static void print_usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -I <genome.idx> -i <reads.fastq> -o <output.sam> -k <max-errors>"
           " [--memory]\n"
        << "\n"
        << "Required arguments:\n"
        << "  -I <path>   Pre-built genome index file\n"
        << "  -i <path>   Input FASTQ reads file\n"
        << "  -o <path>   Output SAM alignment file\n"
        << "  -k <int>    Maximum number of errors (mismatches/edits)\n"
        << "\n"
        << "Options:\n"
        << "  --memory    Load index lazily (lower RAM, slower first access)\n"
        << "  -h, --help  Show this help message\n";
}

// ---- argument parsing (getopt_long) -------------------------------------- //

enum LongOpt : int {
    OPT_MEMORY = 256
};

static MapperConfig parse_args(int argc, char* argv[]) {
    MapperConfig cfg;

    static const char* short_opts = "I:i:o:k:h";

    static struct option long_opts[] = {
        {"memory", no_argument, nullptr, OPT_MEMORY},
        {"help",   no_argument, nullptr, 'h'},
        {nullptr,  0,           nullptr,  0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'I': cfg.idx_path     = optarg; break;
            case 'i': cfg.reads_path   = optarg; break;
            case 'o': cfg.output_path  = optarg; break;
            case 'k': {
                char* endptr = nullptr;
                long val = std::strtol(optarg, &endptr, 10);
                if (endptr == optarg || *endptr != '\0' || val < 0) {
                    std::cerr << "Error: -k requires a non-negative integer.\n";
                    std::exit(EXIT_FAILURE);
                }
                cfg.max_errors = static_cast<int>(val);
                break;
            }
            case OPT_MEMORY:
                cfg.load_mode = LoadMode::MEMORY;
                break;
            case 'h':
                print_usage(argv[0]);
                std::exit(EXIT_SUCCESS);
            default:
                print_usage(argv[0]);
                std::exit(EXIT_FAILURE);
        }
    }

    if (cfg.idx_path.empty() || cfg.reads_path.empty() ||
        cfg.output_path.empty()) {
        std::cerr << "Error: -I, -i, -o, and -k are all required.\n\n";
        print_usage(argv[0]);
        std::exit(EXIT_FAILURE);
    }

    return cfg;
}

struct CandidatePool {
    static constexpr size_t MAX_CANDIDATES = 200;
    std::vector<uint32_t> starts;
    bool anchor_set = false;
    uint32_t anchor = 0;

    CandidatePool() { starts.reserve(MAX_CANDIDATES); }

    void add(uint32_t rs) {
        for (uint32_t x : starts) {
            if (x == rs) return;
        }

        if (!anchor_set) {
            anchor = rs;
            anchor_set = true;
        }

        if (starts.size() < MAX_CANDIDATES) {
            starts.push_back(rs);
            return;
        }

        const auto dist = [](uint32_t a, uint32_t b) -> uint64_t {
            return (a > b) ? static_cast<uint64_t>(a - b)
                           : static_cast<uint64_t>(b - a);
        };

        size_t worst_idx = 0;
        uint64_t worst_dist = dist(starts[0], anchor);
        for (size_t i = 1; i < starts.size(); ++i) {
            const uint64_t d = dist(starts[i], anchor);
            if (d > worst_dist) {
                worst_dist = d;
                worst_idx = i;
            }
        }

        const uint64_t new_dist = dist(rs, anchor);
        if (new_dist < worst_dist)
            starts[worst_idx] = rs;
    }
};

struct Hit {
    int edit = 0;
    uint32_t pos = 0; // 1-based
    std::string cigar;
};

static std::vector<uint8_t> pack_read(const std::string& seq) {
    const size_t bytes = (seq.size() + 3) / 4;
    std::vector<uint8_t> out(bytes, 0);
    for (size_t i = 0; i < seq.size(); ++i)
        out[i / 4] |= static_cast<uint8_t>(dna::ENC_2BIT.v[static_cast<unsigned char>(seq[i])] << ((i % 4) * 2));
    return out;
}

static void push_best_two(const Hit& h, bool& has1, Hit& b1, bool& has2, Hit& b2) {
    if (!has1 || h.edit < b1.edit) {
        if (has1) {
            b2 = b1;
            has2 = true;
        }
        b1 = h;
        has1 = true;
        return;
    }

    if (!has2 || h.edit < b2.edit) {
        b2 = h;
        has2 = true;
    }
}

// ---- main ---------------------------------------------------------------- //

int main(int argc, char* argv[]) {
    MapperConfig cfg = parse_args(argc, argv);

    std::cerr << "[mapper] Index      : " << cfg.idx_path    << "\n"
              << "[mapper] Reads      : " << cfg.reads_path  << "\n"
              << "[mapper] Output     : " << cfg.output_path << "\n"
              << "[mapper] Max errors : " << cfg.max_errors  << "\n"
              << "[mapper] Load mode  : "
              << (cfg.load_mode == LoadMode::MEMORY ? "MEMORY" : "SPEED") << "\n";

    // Load index.
    auto kmer_idx = std::make_unique<KmerIndex>();
    kmer_idx->load(cfg.idx_path, cfg.load_mode);
    std::cerr << "[mapper] Index loaded.\n";

    const std::string chrom = kmer_idx->chrom_name();
    const uint32_t glen = kmer_idx->genome_size();
    if (chrom.empty()) {
        std::cerr << "Error: index does not contain chromosome name metadata.\n";
        return EXIT_FAILURE;
    }

    // Open SAM output.
    SamWriter sam(cfg.output_path);
    sam.write_header(chrom, glen);

    // Map reads.
    FastqReader reader(cfg.reads_path);
    Read r;
    uint64_t total = 0, mapped_count = 0;
    constexpr uint64_t REPORT_EVERY = 100000; // print status every 100k reads
    auto t0 = std::chrono::steady_clock::now();

    while (reader.next(r)) {
        ++total;

        // Uppercase the sequence in place.
        for (char& c : r.seq) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

        const size_t   k        = kmer_idx->k_val();
        const uint32_t read_len = static_cast<uint32_t>(r.seq.size());
        const std::string& read_seq = r.seq;
        const std::string& read_qual = r.qual;

        // --- Pack the read ---
        const auto packed = pack_read(read_seq);

        // --- Handle reads shorter than k: always unmapped ---
        if (read_len < static_cast<uint32_t>(k)) {
            continue;
        }

        // Seed phase: non-overlapping seeds of size index k.
        CandidatePool candidates;
        std::string seed;
        seed.reserve(k);
        for (size_t off = 0; off + k <= read_seq.size(); off += k) {
            seed.assign(read_seq, off, k);
            for (uint32_t gpos : kmer_idx->query(seed)) {
                if (gpos < static_cast<uint32_t>(off)) continue;
                const uint32_t rs = gpos - static_cast<uint32_t>(off);
                if (rs + read_len <= glen)
                    candidates.add(rs);
            }
        }

        bool has1 = false;
        bool has2 = false;
        Hit b1{};
        Hit b2{};

        for (uint32_t rs : candidates.starts) {
            if (cfg.max_errors == 0) {
                // k=0 hot path: avoid substring/DP and use packed SIMD compare.
                const int mm = kmer_idx->verify_match(rs, packed.data(), read_len);
                if (mm == 0) {
                    Hit h;
                    h.edit = 0;
                    h.pos = rs + 1;
                    h.cigar = std::to_string(read_len) + "M";
                    push_best_two(h, has1, b1, has2, b2);
                }
                continue;
            }

            const uint32_t ext_len = std::min<uint32_t>(glen - rs, read_len + static_cast<uint32_t>(cfg.max_errors));
            const std::string text = kmer_idx->genome_substr(rs, ext_len);
            const AlignResult ar = align_banded(read_seq, text, cfg.max_errors);
            if (ar.edit_dist <= cfg.max_errors && !ar.cigar.empty()) {
                Hit h;
                h.edit = ar.edit_dist;
                h.pos = rs + 1;
                h.cigar = ar.cigar;
                push_best_two(h, has1, b1, has2, b2);
                if (has2 && b1.edit == 0 && b2.edit == 0)
                    break;
            }
        }

        if (has1) {
            ++mapped_count;
            if (has2) {
                sam.write_alignment(r.name, chrom, b1.pos, b1.cigar,
                                    read_seq, read_qual,
                                    chrom, b2.pos, b2.cigar);
            } else {
                sam.write_alignment(r.name, chrom, b1.pos, b1.cigar,
                                    read_seq, read_qual);
            }
        }

        if (total % REPORT_EVERY == 0) {
            const auto now = std::chrono::steady_clock::now();
            const double elapsed_s = std::chrono::duration<double>(now - t0).count();
            const double reads_per_s = (elapsed_s > 0.0)
                ? static_cast<double>(total) / elapsed_s
                : 0.0;

            std::cerr << "[mapper] Progress: " << total << " reads, "
                      << mapped_count << " mapped ("
                      << (total ? (100ULL * mapped_count / total) : 0) << "%), "
                      << static_cast<uint64_t>(reads_per_s) << " reads/s, "
                      << static_cast<uint64_t>(elapsed_s) << " s elapsed\n";
        }
    }

    const auto t1 = std::chrono::steady_clock::now();
    const double total_elapsed_s = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "[mapper] Processed " << total   << " reads, "
              << mapped_count << " mapped ("
              << (total ? (100ULL * mapped_count / total) : 0) << "%), "
              << static_cast<uint64_t>(total_elapsed_s) << " s total.\n";

    return EXIT_SUCCESS;
}


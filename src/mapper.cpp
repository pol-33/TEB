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

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_set>
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

// ---- 2-bit read packing -------------------------------------------------- //

static uint8_t base_enc(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return 0;
    }
}

static std::vector<uint8_t> pack_read(const std::string& seq) {
    const size_t bytes = (seq.size() + 3) / 4;
    std::vector<uint8_t> out(bytes, 0);
    for (size_t i = 0; i < seq.size(); ++i)
        out[i / 4] |= static_cast<uint8_t>(base_enc(seq[i]) << ((i % 4) * 2));
    return out;
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
    auto index_base = GenomeIndex::create(IndexAlgo::KMER, OptimizeMode::SPEED);
    auto* kmer_idx  = dynamic_cast<KmerIndex*>(index_base.get());
    if (!kmer_idx) {
        std::cerr << "Error: failed to create KmerIndex.\n";
        return EXIT_FAILURE;
    }
    kmer_idx->load(cfg.idx_path, cfg.load_mode);
    std::cerr << "[mapper] Index loaded.\n";

    // Open SAM output.
    // The genome size and ref name are encoded in the index; we use generic
    // values for the header (SAM still validates; aligners like BWA do the same
    // when ref info isn't readily available).
    SamWriter sam(cfg.output_path);
    sam.write_header("ref", kmer_idx->genome_size());

    // Map reads.
    FastqReader reader(cfg.reads_path);
    Read r;
    uint64_t total = 0, mapped_count = 0;

    while (reader.next(r)) {
        ++total;

        // Uppercase the sequence in place.
        for (char& c : r.seq) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

        const size_t   k        = kmer_idx->k_val();
        const uint32_t read_len = static_cast<uint32_t>(r.seq.size());

        // --- Pack the read ---
        const auto packed = pack_read(r.seq);

        // --- Handle reads shorter than k: always unmapped ---
        if (read_len < static_cast<uint32_t>(k)) {
            sam.write_alignment(r.name, 0, r.seq, 0, false);
            continue;
        }

        // --- Pigeonhole multi-k-mer query ---
        // With ≤ max_errors mismatches in a read, querying k-mers at intervals
        // of  step = floor(k / (max_errors + 1))  guarantees at least one of
        // the queried k-mers is error-free and will be found by the index.
        const size_t step = std::max(size_t{1},
                                     k / static_cast<size_t>(cfg.max_errors + 1));

        // Collect unique candidate read-start positions from all k-mer offsets.
        std::unordered_set<uint32_t> candidate_starts;

        auto add_candidates = [&](size_t off) {
            for (uint32_t gpos : kmer_idx->query(r.seq.substr(off, k))) {
                if (gpos < static_cast<uint32_t>(off)) continue;
                const uint32_t rs = gpos - static_cast<uint32_t>(off);
                if (rs + read_len <= kmer_idx->genome_size())
                    candidate_starts.insert(rs);
            }
        };

        for (size_t off = 0; off + k <= r.seq.size(); off += step)
            add_candidates(off);

        // Always include the last k-mer so the tail of the read is covered.
        const size_t last_off = r.seq.size() - k;
        if (last_off % step != 0)
            add_candidates(last_off);

        // --- Find best-scoring (fewest mismatches) candidate ---
        int      best_mm  = cfg.max_errors + 1;
        uint32_t best_pos = 0;
        bool     found    = false;

        for (uint32_t rs : candidate_starts) {
            const int mm = kmer_idx->verify_match(rs, packed.data(), read_len);
            if (mm <= cfg.max_errors && mm < best_mm) {
                best_mm  = mm;
                best_pos = rs;
                found    = true;
                if (mm == 0) break;  // can't do better
            }
        }

        if (found) {
            ++mapped_count;
            sam.write_alignment(r.name, best_pos, r.seq, best_mm, true);
        } else {
            sam.write_alignment(r.name, 0, r.seq, 0, false);
        }
    }

    std::cerr << "[mapper] Processed " << total   << " reads, "
              << mapped_count << " mapped ("
              << (total ? (100ULL * mapped_count / total) : 0) << "%).\n";

    return EXIT_SUCCESS;
}


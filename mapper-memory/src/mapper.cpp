#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "alignment.hpp"
#include "fasta_reader.hpp"
#include "fastq_reader.hpp"
#include "fm_index.hpp"
#include "fm_search.hpp"
#include "memory_stats.hpp"
#include "nucleotide.hpp"
#include "simd_dispatch.hpp"

// Batch size for parallel processing - tune for cache efficiency
constexpr std::size_t BATCH_SIZE = 1024;

namespace {

struct Config {
    std::string index_path;
    std::string reads_path;
    std::string output_path;
    std::string reference_path;  // Optional: for streaming reference access
    int max_errors = 0;
    int num_threads = 0;  // 0 = auto-detect
};

struct VerifiedAlignment {
    mapper_memory::Alignment alignment;
    std::size_t chrom_index = 0;
    uint64_t text_pos = 0;
};

// Thread-local workspace to avoid heap allocations per read
struct ThreadWorkspace {
    mapper_memory::AlignmentWorkspace align_workspace;
    std::vector<mapper_memory::SearchResult> search_results;
    std::vector<VerifiedAlignment> best_alignments;
    std::unordered_set<uint64_t> seen_positions;
    std::string ref_buffer;
    std::string normalized;
    std::string reverse_complement;
    
    ThreadWorkspace() {
        search_results.reserve(64);
        best_alignments.reserve(2);
        seen_positions.reserve(1024);
    }
};

// Result for a single read - pre-formatted string to avoid locks during output
struct ReadResult {
    std::string output_line;
};

void print_usage() {
    std::cerr << "Usage: mapper -I genome.idx -i reads.fastq -o output.sam -k <0..3> [-t threads] [-R genome.fa]\n"
              << "\n"
              << "Options:\n"
              << "  -I genome.idx    FM-index file (required)\n"
              << "  -i reads.fastq   Input reads file (required)\n"
              << "  -o output.sam    Output file (required)\n"
              << "  -k <0..3>        Maximum edit distance (required)\n"
              << "  -t threads       Number of threads (default: auto-detect)\n"
              << "  -R genome.fa     Reference FASTA (required if index lacks packed genome)\n";
}

Config parse_args(int argc, char* argv[]) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-I" && i + 1 < argc) {
            config.index_path = argv[++i];
        } else if (arg == "-i" && i + 1 < argc) {
            config.reads_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            config.output_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            config.max_errors = std::stoi(argv[++i]);
        } else if (arg == "-t" && i + 1 < argc) {
            config.num_threads = std::stoi(argv[++i]);
        } else if (arg == "-R" && i + 1 < argc) {
            config.reference_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            print_usage();
            std::exit(0);
        } else {
            print_usage();
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (config.index_path.empty() || config.reads_path.empty() || config.output_path.empty() ||
        config.max_errors < 0 || config.max_errors > 3) {
        print_usage();
        throw std::runtime_error("missing required arguments or invalid max_errors");
    }
    return config;
}

bool better_alignment(const VerifiedAlignment& lhs, const VerifiedAlignment& rhs) {
    if (lhs.alignment.edit_dist != rhs.alignment.edit_dist) {
        return lhs.alignment.edit_dist < rhs.alignment.edit_dist;
    }
    if (lhs.alignment.chrom != rhs.alignment.chrom) {
        return lhs.alignment.chrom < rhs.alignment.chrom;
    }
    if (lhs.alignment.ref_pos != rhs.alignment.ref_pos) {
        return lhs.alignment.ref_pos < rhs.alignment.ref_pos;
    }
    return lhs.alignment.cigar < rhs.alignment.cigar;
}

bool same_alignment_bucket(const VerifiedAlignment& lhs,
                           const VerifiedAlignment& rhs,
                           std::size_t read_length) {
    if (lhs.alignment.chrom != rhs.alignment.chrom) {
        return false;
    }
    const uint64_t lhs_pos = lhs.alignment.ref_pos;
    const uint64_t rhs_pos = rhs.alignment.ref_pos;
    const uint64_t delta = (lhs_pos > rhs_pos) ? (lhs_pos - rhs_pos) : (rhs_pos - lhs_pos);
    return delta <= static_cast<uint64_t>(read_length);
}

void maybe_push_alignment(std::vector<VerifiedAlignment>& best,
                          const VerifiedAlignment& candidate,
                          std::size_t read_length) {
    for (VerifiedAlignment& existing : best) {
        if ((existing.chrom_index == candidate.chrom_index && existing.text_pos == candidate.text_pos) ||
            same_alignment_bucket(existing, candidate, read_length)) {
            if (better_alignment(candidate, existing)) {
                existing = candidate;
                std::sort(best.begin(), best.end(), better_alignment);
            }
            return;
        }
    }
    best.push_back(candidate);
    std::sort(best.begin(), best.end(), better_alignment);
    if (best.size() > 2U) {
        best.resize(2U);
    }
}

// Use fast lookup-based complement
std::string reverse_complement_sequence(const std::string& normalized) {
    return mapper_memory::reverse_complement_fast(normalized);
}

mapper_memory::Alignment verify_candidate(const mapper_memory::FMIndexView::ChromosomeView& chrom,
                                          const std::string& normalized_read,
                                          uint64_t text_pos,
                                          int max_errors,
                                          mapper_memory::AlignmentWorkspace& workspace,
                                          std::string& ref_buffer,
                                          mapper_memory::IndexedFasta* fasta) {
    mapper_memory::Alignment best;
    best.edit_dist = std::numeric_limits<int>::max();

    const uint64_t min_ref_len = normalized_read.size() > static_cast<std::size_t>(max_errors)
        ? normalized_read.size() - static_cast<std::size_t>(max_errors)
        : 1U;
    const uint64_t max_ref_len = std::min<uint64_t>(chrom.length - text_pos,
                                                    normalized_read.size() + static_cast<std::size_t>(max_errors));
    if (min_ref_len > max_ref_len) {
        return best;
    }

    ref_buffer.reserve(normalized_read.size() + static_cast<std::size_t>(max_errors));
    for (uint64_t ref_len = min_ref_len; ref_len <= max_ref_len; ++ref_len) {
        // Use either embedded genome or streaming FASTA
        if (chrom.can_extract_reference()) {
            chrom.extract_reference(text_pos, ref_len, ref_buffer);
        } else if (fasta != nullptr) {
            const std::size_t chrom_idx = fasta->find_chromosome(std::string(chrom.name));
            if (chrom_idx != SIZE_MAX) {
                fasta->extract(chrom_idx, text_pos, ref_len, ref_buffer);
            } else {
                continue;
            }
        } else {
            // No way to get reference - skip
            continue;
        }
        
        mapper_memory::Alignment current = mapper_memory::band_align(normalized_read, ref_buffer, max_errors, workspace);
        if (current.edit_dist <= max_errors && current.edit_dist < best.edit_dist) {
            current.chrom = std::string(chrom.name);
            current.ref_pos = text_pos + 1ULL;
            best = current;
        }
    }
    return best;
}

void evaluate_interval(const mapper_memory::FMIndexView::ChromosomeView& chrom,
                       std::size_t chrom_index,
                       const mapper_memory::SearchResult& result,
                       const std::string& normalized_read,
                       int max_errors,
                       mapper_memory::AlignmentWorkspace& workspace,
                       std::string& ref_buffer,
                       std::unordered_set<uint64_t>& seen_positions,
                       std::vector<VerifiedAlignment>& best_alignments,
                       mapper_memory::IndexedFasta* fasta) {
    const uint64_t max_row = std::min<uint64_t>(result.sa_hi, result.sa_lo + 4ULL);
    for (uint64_t row = result.sa_lo; row < max_row; ++row) {
        const uint64_t text_pos = chrom.locate(row);
        if (text_pos >= chrom.length) {
            continue;
        }
        if (!seen_positions.insert((static_cast<uint64_t>(chrom_index) << 32U) ^ text_pos).second) {
            continue;
        }
        if (!chrom.same_chromosome_window(text_pos, normalized_read.size())) {
            continue;
        }
        mapper_memory::Alignment alignment =
            verify_candidate(chrom, normalized_read, text_pos, max_errors, workspace, ref_buffer, fasta);
        if (alignment.edit_dist <= max_errors) {
            maybe_push_alignment(best_alignments, VerifiedAlignment{alignment, chrom_index, text_pos}, normalized_read.size());
        }
    }
}

void search_orientation(const mapper_memory::FMIndexView& index,
                        const std::string& normalized_read,
                        int max_errors,
                        mapper_memory::AlignmentWorkspace& workspace,
                        std::string& ref_buffer,
                        std::unordered_set<uint64_t>& seen_positions,
                        std::vector<mapper_memory::SearchResult>& search_results,
                        std::vector<VerifiedAlignment>& best_alignments,
                        mapper_memory::IndexedFasta* fasta) {
    for (std::size_t chrom_index = 0; chrom_index < index.chromosome_count(); ++chrom_index) {
        const auto& chrom = index.chromosome(chrom_index);
        search_results.clear();

        const mapper_memory::SAInterval exact = mapper_memory::exact_search(chrom, normalized_read);
        if (exact.lo < exact.hi) {
            search_results.push_back(mapper_memory::SearchResult{exact.lo, exact.hi, 0});
        } else if (max_errors > 0) {
            mapper_memory::inexact_search(chrom, normalized_read, max_errors, search_results, 32U);
        }
        if (search_results.empty()) {
            continue;
        }

        const int best_interval_edit = search_results.front().edit_dist;
        for (const mapper_memory::SearchResult& result : search_results) {
            if (result.edit_dist > best_interval_edit) {
                break;
            }
            evaluate_interval(chrom, chrom_index, result, normalized_read, max_errors, workspace,
                              ref_buffer, seen_positions, best_alignments, fasta);
            if (best_alignments.size() == 2U && best_alignments.back().alignment.edit_dist == 0) {
                break;
            }
        }
    }
}

// Thread-safe version that returns a string
std::string format_record(const mapper_memory::Read& read,
                          const std::vector<VerifiedAlignment>& alignments) {
    std::string result;
    result.reserve(512);
    
    if (alignments.empty()) {
        result += read.name;
        result += " * 0 * ";
        result += read.seq;
        result += ' ';
        result += read.qual;
        result += '\n';
        return result;
    }
    
    const mapper_memory::Alignment& primary = alignments.front().alignment;
    result += read.name;
    result += ' ';
    result += primary.chrom;
    result += ' ';
    result += std::to_string(primary.ref_pos);
    result += ' ';
    result += primary.cigar;
    result += ' ';
    result += read.seq;
    result += ' ';
    result += read.qual;
    
    if (alignments.size() > 1U) {
        const mapper_memory::Alignment& alt = alignments[1].alignment;
        result += " ALT:";
        result += alt.chrom;
        result += ',';
        result += std::to_string(alt.ref_pos);
        result += ',';
        result += alt.cigar;
    }
    result += '\n';
    return result;
}

void report_peak_rss() {
    const mapper_memory::MemoryStats stats = mapper_memory::read_memory_stats();
    std::cerr << "[mapper] current RSS: " << mapper_memory::bytes_to_mebibytes(stats.current_rss_bytes)
              << " MiB, peak RSS: " << mapper_memory::bytes_to_mebibytes(stats.peak_rss_bytes) << " MiB\n";
}

// Process a single read using thread-local workspace
std::string process_read(const mapper_memory::FMIndexView& index,
                         const mapper_memory::Read& read,
                         int max_errors,
                         ThreadWorkspace& ws,
                         mapper_memory::IndexedFasta* fasta) {
    ws.normalized = mapper_memory::normalize_sequence(read.seq);
    ws.best_alignments.clear();
    ws.seen_positions.clear();

    search_orientation(index, ws.normalized, max_errors, ws.align_workspace, ws.ref_buffer,
                       ws.seen_positions, ws.search_results, ws.best_alignments, fasta);

    ws.reverse_complement = reverse_complement_sequence(ws.normalized);
    if (ws.reverse_complement != ws.normalized || ws.best_alignments.empty()) {
        search_orientation(index, ws.reverse_complement, max_errors, ws.align_workspace, ws.ref_buffer,
                           ws.seen_positions, ws.search_results, ws.best_alignments, fasta);
    }

    return format_record(read, ws.best_alignments);
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        std::ios::sync_with_stdio(false);
        const Config config = parse_args(argc, argv);

        // Determine number of threads
        int num_threads = config.num_threads;
        if (num_threads <= 0) {
            num_threads = static_cast<int>(std::thread::hardware_concurrency());
            if (num_threads <= 0) num_threads = 1;
        }

#ifdef _OPENMP
        omp_set_num_threads(num_threads);
        std::cerr << "[mapper] using " << num_threads << " threads (OpenMP)\n";
#else
        if (num_threads > 1) {
            std::cerr << "[mapper] warning: compiled without OpenMP, using 1 thread\n";
        }
        num_threads = 1;
#endif

        std::cerr << "[mapper] loading FM-index from " << config.index_path << "\n";
        mapper_memory::FMIndexView index(config.index_path);
        std::cerr << "[mapper] index loaded, " << index.chromosome_count() << " chromosomes\n";
        const mapper_memory::simd::DispatchInfo simd_dispatch = mapper_memory::simd::resolved_dispatch();
        std::cerr << "[mapper] occ backend: " << mapper_memory::simd::backend_name(simd_dispatch.backend) << "\n";
        
        // Check if we need streaming FASTA
        std::unique_ptr<mapper_memory::IndexedFasta> streaming_fasta;
        const bool index_has_genome = index.has_genome() && 
            (index.chromosome_count() == 0 || index.chromosome(0).can_extract_reference());
        
        if (!index_has_genome) {
            if (config.reference_path.empty()) {
                throw std::runtime_error("index does not contain packed genome; please provide -R genome.fa");
            }
            std::cerr << "[mapper] index lacks packed genome, opening streaming FASTA from " << config.reference_path << "\n";
            streaming_fasta = std::make_unique<mapper_memory::IndexedFasta>(config.reference_path);
            std::cerr << "[mapper] streaming FASTA opened, " << streaming_fasta->chromosome_count() << " chromosomes\n";
        }
        
        std::cerr << "[mapper] opening reads from " << config.reads_path << "\n";
        mapper_memory::FastqReader reader(config.reads_path);
        std::ofstream out(config.output_path);
        if (!out) {
            throw std::runtime_error("failed to open output file: " + config.output_path);
        }
        std::vector<char> output_buffer(4 * 1024 * 1024);
        out.rdbuf()->pubsetbuf(output_buffer.data(), static_cast<std::streamsize>(output_buffer.size()));

        std::cerr << "[mapper] starting read processing with k=" << config.max_errors << "\n";

        // Create thread-local workspaces (minimal additional memory)
        std::vector<ThreadWorkspace> workspaces(static_cast<std::size_t>(num_threads));

        // Batch storage
        std::vector<mapper_memory::Read> batch_reads;
        std::vector<std::string> batch_results;
        batch_reads.reserve(BATCH_SIZE);
        batch_results.resize(BATCH_SIZE);

        uint64_t processed = 0;
        mapper_memory::Read read;
        
        // Get raw pointer for use in parallel region
        mapper_memory::IndexedFasta* fasta_ptr = streaming_fasta.get();

        while (true) {
            // Read a batch (sequential - fast I/O)
            batch_reads.clear();
            while (batch_reads.size() < BATCH_SIZE && reader.next(read)) {
                batch_reads.push_back(read);
            }
            
            if (batch_reads.empty()) {
                break;
            }

            const std::size_t batch_size = batch_reads.size();
            batch_results.resize(batch_size);

            // Process batch in parallel
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
#endif
            for (std::size_t i = 0; i < batch_size; ++i) {
#ifdef _OPENMP
                const int tid = omp_get_thread_num();
#else
                const int tid = 0;
#endif
                batch_results[i] = process_read(index, batch_reads[i], config.max_errors, 
                                                 workspaces[static_cast<std::size_t>(tid)], fasta_ptr);
            }

            // Write results in order (sequential - maintains ordering)
            for (std::size_t i = 0; i < batch_size; ++i) {
                out << batch_results[i];
            }

            processed += batch_size;

            // Progress reporting - show after first batch, then at milestones
            if (processed <= BATCH_SIZE || processed % 10000ULL < BATCH_SIZE) {
                std::cerr << "[mapper] processed " << processed << " reads\n";
            }
        }

        std::cerr << "[mapper] processed " << processed << " reads total\n";
        report_peak_rss();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mapper error: " << e.what() << '\n';
        return 1;
    }
}

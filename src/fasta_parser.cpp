#include "fasta_parser.hpp"
#include <sys/mman.h>   // mmap / munmap / madvise
#include <sys/stat.h>   // fstat
#include <fcntl.h>      // open
#include <unistd.h>     // close

#define KMER_INDEX_SIZE 10

using namespace std;

// mmap-based single-pass FASTA parser.
// Uses memchr to find newlines (~32 bytes/cycle via SIMD libc), count_gc_bulk
// directly on the mapped buffer (zero copies), and a 1 MB write buffer for output.
static void parseFastaFile(const string& infile, const string& outfile, GlobalStats& stats, const unsigned int kmer_length, bool per_seq_stats) {
    // Memory-map the input file
    int fd = open(infile.c_str(), O_RDONLY);
    if (fd < 0) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    }

    // Get file size and map it into memory
    struct stat sb;
    if (fstat(fd, &sb) < 0) {
        cerr << "Error: fstat failed for " << infile << endl;
        close(fd); return;
    }
    size_t filesize = (size_t)sb.st_size;
    
    const char* data = static_cast<const char*>(
        mmap(nullptr, filesize, PROT_READ, MAP_PRIVATE, fd, 0));

    close(fd);  // fd no longer needed after mmap

    if (data == MAP_FAILED) {
        cerr << "Error: mmap failed for " << infile << endl;
        return;
    }
    // Tell the kernel our access pattern.
    // Tell the kernel we will read this linearly, so we can read ahead
    madvise(const_cast<char*>(data), filesize, MADV_SEQUENTIAL);


    static char out_buf[IO_BUFFER_SIZE];
    ofstream out;
    if (!outfile.empty()) {
        out.rdbuf()->pubsetbuf(out_buf, sizeof(out_buf));
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
            // Free the mapped memory before returning
            munmap(const_cast<char*>(data), filesize);
            return;
        }
    }

    kmer_table_t kmer_indxs;

    const char* p   = data;
    const char* end = data + filesize;

    string current_header;
    SequenceStats cur_seq = {0, 0, 0.0};

    // Finalise current sequence: only populate per_sequence map when requested.
    auto finalizeSequence = [&]() {
        if (!current_header.empty() && per_seq_stats) {
            cur_seq.gc_content = (cur_seq.length > 0)
                ? (double)cur_seq.gc_count / cur_seq.length * 100.0 : 0.0;
            stats.per_sequence[current_header] = cur_seq;
        }
    };

    while (p < end) {
        // Find next newline, starts at p and ends at nl (if nl is null, end of file reached)
        const char* nl = static_cast<const char*>(memchr(p, '\n', (size_t)(end - p)));
        const char* line_end = nl ? nl : end;
        size_t line_len = (size_t)(line_end - p);

        // Strip trailing \r (Windows / mixed line endings)
        if (line_len && p[line_len - 1] == '\r') --line_len;

        if (line_len == 0) { p = nl ? nl + 1 : end; continue; }  // skip blank lines

        if (*p == '>') {
            finalizeSequence();
            if (!current_header.empty()) out.put('\n');  // Finish previous record with a newline
            current_header.assign(p + 1, line_len - 1);  // text after '>'
            cur_seq = {0, 0, 0.0};
            ++stats.num_sequences;
            if (out.is_open()) {
                out.write(p, (streamsize)line_len);  // includes '>'
                out.put('\n');
            }
        } else {
            // Sequence line: bulk GC directly on mapped buffer
            long long gc = count_gc_bulk(p, line_len);
            stats.total_length   += (long long)line_len;
            stats.total_gc_count += gc;
            cur_seq.length       += (long long)line_len;
            cur_seq.gc_count     += gc;
            if (out.is_open()) {
                // Concatenate sequence lines until the next header
                out.write(p, (streamsize)line_len);
            }
            if (kmer_length > 0)
                update_kmer_table(string(p, line_len), kmer_indxs, kmer_length);
        }

        p = nl ? nl + 1 : end;
    }
    finalizeSequence();  // store the last sequence

    stats.overall_gc_content = (stats.total_length > 0)
        ? (double)stats.total_gc_count / stats.total_length * 100.0 : 0.0;

    // Free the mapped memory
    munmap(const_cast<char*>(data), filesize);

    #ifdef DEBUG
    print_kmer_table(kmer_indxs);
    #endif
}

// ─────────────────────────────────────────────────────────────────────────────
// parseFastaFileStream  —  low-memory mode (-m flag)
//
// Uses read() with a 64 MB rolling buffer instead of mmap so that only the
// buffer lives in the process's RSS, not the entire file.
// ─────────────────────────────────────────────────────────────────────────────
static void parseFastaFileStream(const string& infile, const string& outfile,
                                  GlobalStats& stats, const unsigned int kmer_length, bool per_seq_stats, size_t stream_buf) {
    int fd = open(infile.c_str(), O_RDONLY);
    if (fd < 0) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    }

    vector<char> buf(stream_buf);

    static char out_buf[IO_BUFFER_SIZE];
    ofstream out;
    if (!outfile.empty()) {
        out.rdbuf()->pubsetbuf(out_buf, sizeof(out_buf));
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Error: The file " << outfile << " could not be opened/created." << endl;
            close(fd); return;
        }
    }

    kmer_table_t kmer_indxs;
    string current_header;
    SequenceStats cur_seq = {0, 0, 0.0};

    auto finalizeSequence = [&]() {
        if (!current_header.empty() && per_seq_stats) {
            cur_seq.gc_content = (cur_seq.length > 0)
                ? (double)cur_seq.gc_count / cur_seq.length * 100.0 : 0.0;
            stats.per_sequence[current_header] = cur_seq;
        }
    };

    size_t filled = 0;
    bool   eof    = false;

    while (!eof || filled > 0) {
        if (!eof) {
            ssize_t n = read(fd, buf.data() + filled, stream_buf - filled);
            if (n <= 0) { eof = true; }
            else        { filled += (size_t)n; }
        }

        const char* base = buf.data();
        const char* p    = base;
        const char* end  = base + filled;

        while (p < end) {
            const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
            if (!nl) {
                if (eof) nl = end;  // last line without trailing newline
                else break;         // need more data
            }
            const char* line_end = nl;
            size_t line_len = (size_t)(line_end - p);
            if (line_len && p[line_len - 1] == '\r') --line_len;

            if (line_len == 0) { p = nl < end ? nl + 1 : end; continue; }

            if (*p == '>') {
                finalizeSequence();
                if (!current_header.empty()) out.put('\n');  // close previous sequence
                current_header.assign(p + 1, line_len - 1);
                cur_seq = {0, 0, 0.0};
                ++stats.num_sequences;
                if (out.is_open()) {
                    out.write(p, (streamsize)line_len);
                    out.put('\n');
                }
            } else {
                long long gc = count_gc_bulk(p, line_len);
                stats.total_length   += (long long)line_len;
                stats.total_gc_count += gc;
                cur_seq.length       += (long long)line_len;
                cur_seq.gc_count     += gc;
                if (out.is_open()) {
                    out.write(p, (streamsize)line_len);
                }
                if (kmer_length > 0)
                    update_kmer_table(string(p, line_len), kmer_indxs, kmer_length);
            }

            p = nl < end ? nl + 1 : end;
        }

        // Compact: shift unprocessed tail to start of buffer
        size_t leftover = (size_t)(end - p);
        if (leftover > 0 && p != base)
            memmove(buf.data(), p, leftover);
        filled = leftover;

        if (!eof && filled == stream_buf)
            throw runtime_error("[parseFastaFileStream] line exceeds buffer size");
    }

    finalizeSequence();

    stats.overall_gc_content = (stats.total_length > 0)
        ? (double)stats.total_gc_count / stats.total_length * 100.0 : 0.0;

    close(fd);

#ifdef DEBUG
    print_kmer_table(kmer_indxs);
#endif
}

// Print the statistics
void printStatistics(const GlobalStats& stats) {
    cout << "============================================" << endl;
    cout << "       SUMMARY GENOMIC ANALYSIS" << endl;
    cout << "============================================" << endl;
    cout << "Number of processed sequences: " << stats.num_sequences << endl;

    // Per-sequence block is printed only when the map was populated (per_seq_stats=true).
    if (!stats.per_sequence.empty()) {
        cout << "--------------------------------------------" << endl;
        for (const auto& [id, seq] : stats.per_sequence) {
            cout << "ID: " << id << endl;
            cout << "  > Length: " << seq.length << " bp" << endl;
            cout << "  > GC content: " << fixed << setprecision(2) << seq.gc_content << "%" << endl;
            cout << endl;
            cout << "--------------------------------------------" << endl;
        }
    }

    cout << "GLOBAL SUMMARY:" << endl;
    cout << "  > Total length: " << stats.total_length << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "============================================" << endl;
}

int fasta_parser(const string& input_file, const string& output_file, const unsigned int kmer_length, bool per_seq_stats, bool low_mem, size_t stream_buf) {
    GlobalStats stats;
    stats.num_sequences = 0;
    stats.total_length = 0;
    stats.total_gc_count = 0;

    if (low_mem)
        parseFastaFileStream(input_file, output_file, stats, kmer_length, per_seq_stats, stream_buf);
    else
        parseFastaFile(input_file, output_file, stats, kmer_length, per_seq_stats);
    if (stats.num_sequences == 0) {
        cout << "No sequences found." << endl;
        return 0;
    }

    printStatistics(stats);

    return 0;
}


#include "fastq_parser.hpp"
#include "utils.hpp"

#include <climits>
#include <iomanip>
#include <fcntl.h>      // open
#include <unistd.h>     // close

using namespace std;


// parseFastqFileStream by streaming the file
//
// Uses read() with a configurable rolling buffer. Peak RSS is bounded to the
// buffer size regardless of file size.
static void parseFastqFileStream(const string& infile, const string& outfile,
                                  FastqGlobalStats& stats, const int qmin, bool per_seq_stats, size_t stream_buf) {
    int fd = open(infile.c_str(), O_RDONLY);
    if (fd < 0) {
        cerr << "Error: The file " << infile << " could not be opened." << endl;
        return;
    }

    // Ask the kernel for aggressive read-ahead on this fd.
#ifdef __APPLE__
    fcntl(fd, F_RDAHEAD, 1);
#else
    posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

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

    stats.num_sequences  = 0;
    stats.total_length   = 0;
    stats.total_gc_count = 0;
    stats.minimum        = LLONG_MAX;
    stats.maximum        = 0LL;

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
            // Skip blank / non-@ lines
            while (p < end && (unsigned char)*p <= (unsigned char)'\r') ++p;
            if (p >= end) break;

            if (*p != '@') {
                const char* nl = (const char*)memchr(p, '\n', (size_t)(end - p));
                p = nl ? nl + 1 : end;
                continue;
            }

            // Line 1: header
            const char* nl1 = (const char*)memchr(p, '\n', (size_t)(end - p));
            if (!nl1) break;
            // Line 2: sequence
            const char* seq_start = nl1 + 1;
            const char* nl2 = (const char*)memchr(seq_start, '\n', (size_t)(end - seq_start));
            if (!nl2) break;
            // Line 3: '+'
            const char* sep_start = nl2 + 1;
            const char* nl3 = (const char*)memchr(sep_start, '\n', (size_t)(end - sep_start));
            if (!nl3) break;
            // Line 4: quality (may lack trailing newline at EOF)
            const char* qual_start = nl3 + 1;
            const char* nl4 = (const char*)memchr(qual_start, '\n', (size_t)(end - qual_start));
            if (!nl4) {
                if (!eof) break;   // incomplete record — need more data
                nl4 = end;         // last record without trailing newline
            }

            const char* hdr_start = p + 1;
            size_t hdr_len  = (size_t)(nl1 - hdr_start);
            if (hdr_len  && hdr_start[hdr_len - 1]  == '\r') --hdr_len;
            size_t seq_len  = (size_t)(nl2 - seq_start);
            if (seq_len  && seq_start[seq_len - 1]  == '\r') --seq_len;
            size_t qual_len = (size_t)(nl4 - qual_start);
            if (qual_len && qual_start[qual_len - 1] == '\r') --qual_len;

            if (seq_len != qual_len)
                throw runtime_error("[parseFastqFileStream] sequence length != quality length");

            size_t keep = seq_len;
            if (qmin > 0) keep = get_trim_limit(qual_start, qual_len, qmin);

            if (out.is_open()) {
                out.put('@');
                out.write(hdr_start,  (streamsize)hdr_len);  out.put('\n');
                out.write(seq_start,  (streamsize)keep);     out.put('\n');
                out.put('+');                                 out.put('\n');
                out.write(qual_start, (streamsize)keep);     out.put('\n');
            }

            long long gc  = count_gc_bulk(seq_start, keep);
            long long len = (long long)keep;
            stats.total_length   += len;
            stats.total_gc_count += gc;
            if (len < stats.minimum) stats.minimum = len;
            if (len > stats.maximum) stats.maximum = len;
            ++stats.num_sequences;

            if (per_seq_stats) {
                SequenceStats sStats;
                sStats.length     = len;
                sStats.gc_count   = gc;
                sStats.gc_content = (len > 0) ? (double)gc / len * 100.0 : 0.0;
                stats.per_sequence.emplace(string(hdr_start, hdr_len), sStats);
            }

            p = (nl4 < end) ? nl4 + 1 : end;
        }

        // Compact: shift unprocessed tail to start of buffer
        size_t leftover = (size_t)(end - p);
        if (leftover > 0 && p != base)
            memmove(buf.data(), p, leftover);
        filled = leftover;

        if (!eof && filled == stream_buf)
            throw runtime_error("[parseFastqFileStream] record exceeds buffer size");
    }

    close(fd);

    stats.overall_gc_content = (stats.total_length > 0)
        ? (double)stats.total_gc_count / stats.total_length * 100.0 : 0.0;
    stats.avg_read_len = (stats.num_sequences > 0)
        ? stats.total_length / stats.num_sequences : 0;
}

// Print the statistics
void printStatistics(const FastqGlobalStats& stats) {
    cout << "============================================" << endl;
    cout << "       SUMMARY GENOMIC ANALYSIS" << endl;
    cout << "============================================" << endl;
    cout << "Number of processed sequences: " << stats.num_sequences << endl;
    if (!stats.per_sequence.empty()) {
        cout << "--------------------------------------------" << endl;
        for (auto const& [header, sStats] : stats.per_sequence) {
            cout << "ID: " << header << endl;
            cout << "  > Length: " << sStats.length << " bp" << endl;
            cout << "  > GC content: " << fixed << setprecision(2) << sStats.gc_content << "%" << endl;
            cout << endl;
        }
    }
    cout << "--------------------------------------------" << endl;
    cout << "GLOBAL SUMMARY:" << endl;
    cout << "  > Total length: "     << stats.total_length   << " bp" << endl;
    cout << "  > Total GC content: " << fixed << setprecision(2) << stats.overall_gc_content << "%" << endl;
    cout << "  > Minimum read: "     << stats.minimum        << " bp" << endl;
    cout << "  > Maximum read: "     << stats.maximum        << " bp" << endl;
    cout << "  > Avg length read: "  << stats.avg_read_len   << " bp" << endl;
    cout << "============================================" << endl;
}

int fastq_parser(const string& input_file, const string& output_file, const int qmin, bool per_seq_stats, size_t stream_buf) {

    FastqGlobalStats stats;

    parseFastqFileStream(input_file, output_file, stats, qmin, per_seq_stats, stream_buf);

    if (stats.num_sequences == 0) {
        cout << "No sequences found." << endl;
        return 0;
    }

    printStatistics(stats);
    return 0;
}


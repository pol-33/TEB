# mapper-memory Documentation

## 1. Overview

`mapper-memory` is the memory-oriented implementation of the final sequence-mapping project for
the course *Bioinformatics Algorithms and Techniques*. Its goal is to align short reads from a
FASTQ file against a reference genome from a FASTA file, allowing a maximum edit distance
`k` with `0 <= k <= 3`, and to report at most two best alignments for each read.

This mapper is built around a different priority from `mapper-speed`:

- minimize index size,
- minimize peak RSS during mapping,
- remain correct,
- still keep the implementation practical enough to benchmark on large references such as GRCh38.

The design therefore replaces the large direct-address seed index with a per-chromosome
FM-index. This reduces the persistent index footprint dramatically, at the cost of more
expensive query-time navigation.

## 2. Project Interface

The project still follows the required two-program interface.

### Index construction

```bash
./indexer -R genome.fa -I genome.fmidx
```

### Mapping

```bash
./mapper -I genome.fmidx -i reads.fastq -o output.sam -k 1 -R genome.fa
```

The `indexer` executable expects:

- `-R`: reference FASTA path,
- `-I`: output FM-index path,
- `-L`: optional memory optimization level from `0` to `5`.

The `mapper` executable expects:

- `-I`: FM-index path,
- `-i`: input FASTQ path,
- `-o`: output SAM-like path,
- `-k`: maximum edit distance from `0` to `3`,
- `-t`: optional thread count,
- `-R`: reference FASTA path when the index does not embed the packed genome.

## 3. Default and Recommended Mode

The default indexing level is now `5`, which is the most aggressive memory-optimization mode.

This means the default command:

```bash
./indexer -R genome.fa -I genome.fmidx
```

is equivalent to:

```bash
./indexer -R genome.fa -I genome.fmidx -L 5
```

At level `5`:

- Occ checkpoints are sparse,
- SA samples are sparse,
- the packed genome is not embedded in the FM-index,
- the mapper must therefore receive `-R genome.fa`.

So the default memory-contest workflow is:

```bash
./indexer -R genome.fa -I genome.fmidx
./mapper -I genome.fmidx -i reads.fastq -o output.sam -k 1 -R genome.fa
```

If the goal is not minimum memory but rather a speed-oriented baseline, level `0` remains
available:

```bash
./indexer -R genome.fa -I genome.fmidx -L 0
./mapper -I genome.fmidx -i reads.fastq -o output.sam -k 1
```

## 4. General Design Philosophy

The pipeline is FM-index based and chromosome segmented:

1. Read the FASTA chromosome by chromosome.
2. Build one suffix array per chromosome.
3. Convert that suffix array into one chromosome-local packed BWT.
4. Store Occ checkpoints, sampled SA values, and optionally a packed chromosome sequence.
5. Memory-map the serialized FM-index at runtime.
6. Search reads with exact backward search first.
7. Fall back to bounded inexact FM backtracking when exact search fails.
8. Recover candidate genomic positions with `locate()`.
9. Verify candidates with a banded edit-distance DP.
10. Report the best two simplified SAM-like hits.

This structure pushes most persistent memory into compressed, mmap-friendly structures and keeps
runtime workspace modest.

## 5. Sequence Utilities and Core Encodings

The low-level sequence utilities live mainly in:

- `include/nucleotide.hpp`
- `include/memory_stats.hpp`

### 5.1 Nucleotide encoding

`nucleotide.hpp` defines a compact six-symbol alphabet:

- sentinel,
- separator,
- A,
- C,
- G,
- T.

It also provides:

- fast ASCII-based normalization,
- fast base-to-rank conversion,
- reverse-complement helpers,
- packed 2-bit storage helpers.

An important design choice is that ambiguous bases are normalized into the DNA alphabet rather
than being stored in a separate ambiguity structure. This is simpler and smaller than the
`mapper-speed` representation, but it also means this mapper is more aggressively simplified in
how it treats ambiguous reference or read symbols.

### 5.2 Memory accounting

`memory_stats.hpp` provides portable RSS accounting:

- on Linux, by reading `/proc/self/status` and falling back to `getrusage`,
- on macOS, by using Mach task info plus `getrusage`.

These helpers are used by both `indexer` and `mapper` to print current and peak RSS to stderr.

## 6. FASTA and FASTQ Input

The input layer lives in:

- `include/fasta_reader.hpp`
- `src/fasta_reader.cpp`
- `include/fastq_reader.hpp`
- `src/fastq_reader.cpp`

### 6.1 Sequential FASTA reader

`FastaReader` streams one chromosome at a time from a FASTA file.

Its responsibilities are:

- normalize chromosome names,
- strip CRLF if necessary,
- normalize bases to uppercase A/C/G/T,
- expose one `FastaChromosome` at a time.

This is used by `indexer`, which never needs the whole reference as a single concatenated object.

### 6.2 Indexed streaming FASTA

`IndexedFasta` is the more memory-sensitive component.

It builds or loads a `.fai`-style index that stores:

- chromosome name,
- chromosome length,
- byte offset of the first base,
- number of bases per FASTA line,
- number of bytes per FASTA line.

With that metadata, it can extract only the reference window needed for candidate verification.
This is the mechanism that lets levels `1` to `5` omit the packed genome from the FM-index file.

### 6.3 FASTQ reader

`FastqReader` is deliberately small and strict:

- it validates the `@` header line,
- it validates the `+` separator line,
- it checks sequence/quality length consistency,
- it trims header names at the first whitespace.

This keeps mapping logic simple and gives early failures on malformed input.

## 7. Suffix Array and BWT Construction

The suffix-array and BWT logic lives in:

- `include/bwt.hpp`
- `src/bwt.cpp`

### 7.1 SA construction strategy

`src/bwt.cpp` implements a layered suffix-array builder:

- for very small inputs, a naive comparator-based suffix sort,
- for small-to-medium inputs, a doubling algorithm,
- for larger inputs, SA-IS.

This keeps the implementation self-contained while still scaling to chromosome-sized texts.

### 7.2 Why chromosome-local SA instead of whole-genome SA

The most important memory decision in this module is that it builds one FM-index segment per
chromosome instead of a single global FM-index over the whole reference.

This helps because:

- suffix-array construction stays bounded by chromosome size,
- temporary memory spikes are smaller,
- indexing failures are less likely on machines with limited RAM,
- the final serialized format remains naturally segmented.

### 7.3 BWT payload

`build_bwt()` produces:

- the text length,
- the primary index,
- rank counts,
- the packed BWT,
- separator-row bookkeeping.

The packed BWT uses the same 2-bit DNA packing style as the rest of the mapper.

## 8. FM-index Layout and Serialization

The FM-index interface and implementation live in:

- `include/fm_index.hpp`
- `src/fm_index.cpp`

### 8.1 Per-chromosome representation

Each `OwnedChromosomeIndex` stores:

- chromosome name,
- chromosome length,
- text length,
- BWT primary row,
- Occ checkpoint stride,
- SA sample stride,
- cumulative count array,
- packed BWT,
- sampled Occ array,
- sampled SA array,
- optional packed genome.

At runtime, `FMIndexView` maps the whole file into memory and creates lightweight
`ChromosomeView` objects that point directly into the mapped bytes.

### 8.2 Serialized file format

The on-disk format stores:

- a magic/version header,
- global Occ and SA strides,
- chromosome count,
- feature flags,
- one header per chromosome,
- aligned binary payload blocks for packed BWT, Occ, SA, and optional genome.

Version `002` supports optional per-chromosome genome embedding and remains backward compatible
with the earlier `FMCHR001` format.

### 8.3 FM operations

The key runtime FM-index operations are:

- `occ(rank, pos)`,
- `lf(row)`,
- `locate(row)`.

`occ()` combines:

- the nearest stored checkpoint,
- a local count over the packed BWT suffix between that checkpoint and the query position.

`locate()` repeatedly applies LF-mapping until it reaches a sampled SA row.

This is precisely why sparser SA and Occ sampling reduce memory but slow mapping.

## 9. Optimization Levels 0 to 5

The optimization-level definitions are in `include/fm_index.hpp`, and the CLI mapping is handled
in `src/indexer.cpp`.

The levels are:

| Level | Occ stride | SA stride | Store genome | Main effect |
| --- | ---: | ---: | --- | --- |
| 0 | 256 | 32 | yes | largest, fastest |
| 1 | 256 | 32 | no | remove embedded genome |
| 2 | 256 | 128 | no | sparse SA |
| 3 | 512 | 256 | no | sparse SA + sparse Occ |
| 4 | 1024 | 512 | no | very sparse |
| 5 | 2048 | 1024 | no | maximum memory reduction |

Interpretation:

- increasing the Occ stride reduces checkpoint memory but increases the scan work inside `occ()`,
- increasing the SA stride reduces sample memory but increases LF steps in `locate()`,
- disabling embedded genome removes index bytes but forces reference windows to come from `-R`.

### 9.1 Why level 5 is now the default

The module exists for the memory contest. Making level `5` the default aligns the default user
experience with that objective:

- smaller serialized index,
- lower mapping RSS,
- no need to remember `-L 5` explicitly.

The tradeoff is very large runtime cost.

The archived benchmark outputs in:

- `benchmark-dataset-complet/summary.txt`
- `benchmark-dataset-complet-max-level/summary.txt`

show that this tradeoff is real:

- the baseline configuration reached roughly `1.88 GiB` peak RSS on the archived run,
- the max-level configuration dropped to roughly `0.78 GiB`,
- but runtime also became dramatically slower.

So `L5` is the right default for memory-first benchmarking, not for throughput-first runs.

## 10. Exact and Inexact Search

The search logic lives in:

- `include/fm_search.hpp`
- `src/fm_search.cpp`

### 10.1 Exact search

`exact_search()` performs classic backward search over one chromosome-local FM-index.

For each read character from right to left, it updates:

- `lo = C[c] + Occ(c, lo)`
- `hi = C[c] + Occ(c, hi)`

If the interval becomes empty, exact search fails for that chromosome.

### 10.2 Inexact search

`inexact_search()` uses recursive backtracking with lower-bound pruning.

The recursion supports:

- exact match,
- substitution,
- deletion,
- insertion.

To keep the search bounded:

- a lower-bound array estimates how many edit segments are already forced,
- recursion stops once `max_results` is reached,
- the mapper later only evaluates the best edit-distance intervals first.

This is still far more expensive than exact search, which is why FM-index navigation is a major
runtime bottleneck in memory-optimized levels.

## 11. Candidate Verification and Alignment

Candidate verification lives in:

- `include/alignment.hpp`
- `src/alignment.cpp`
- `src/mapper.cpp`

### 11.1 Reference window recovery

Once a search interval is found, the mapper:

- uses `locate()` to recover candidate text positions,
- deduplicates candidate positions,
- checks that the read-length window stays inside the chromosome,
- recovers reference text either from the embedded packed genome or from `IndexedFasta`.

### 11.2 Banded edit-distance DP

The verifier is a reusable banded dynamic program:

- insertions, deletions, and substitutions are supported,
- traceback reconstructs a simplified CIGAR over `M`, `I`, and `D`,
- only windows compatible with `k` are explored.

This verifier is smaller and simpler than the SIMD-heavy verifier in `mapper-speed`. Here the
main optimization target is not DP throughput but limiting how often verification must happen.

### 11.3 Reporting policy

The mapper reports:

- unmapped reads as `* 0 *`,
- up to two alignments,
- one primary alignment plus one `ALT:` alignment when a second candidate survives.

## 12. Mapper Execution Model

The top-level mapper implementation is `src/mapper.cpp`.

Important runtime choices:

- reads are processed in batches of `1024`,
- read input and output are sequential,
- batch work is parallelized with OpenMP when available,
- each thread owns its own `ThreadWorkspace`,
- output order is preserved by writing completed batch results sequentially.

This means the mapper is parallel over reads, not over one read’s search itself.

### 12.1 Why thread-local workspaces matter

Each thread keeps:

- one DP workspace,
- one search-results buffer,
- one best-alignments buffer,
- one seen-position hash set,
- temporary normalized/ref buffers.

This avoids per-read heap churn and keeps peak additional RAM from threading relatively modest
compared with the mmap-backed FM-index itself.

## 13. SIMD Dispatch

The SIMD backend lives in:

- `include/simd_dispatch.hpp`
- `src/simd_dispatch.cpp`
- `src/simd_avx512.cpp`
- `src/simd_selftest.cpp`

### 13.1 What SIMD accelerates

In this mapper, SIMD only accelerates the packed-BWT counting inside `count_packed_range()`.

That means AVX-512 helps only when a substantial fraction of runtime is spent in:

- `occ()`,
- `lf()`,
- `locate()`,
- exact and inexact FM navigation.

It does not accelerate the full mapper uniformly.

### 13.2 Scalar backend

The scalar backend uses:

- a precomputed per-byte lookup table,
- partial-byte handling at the range boundaries,
- explicit correction for the sentinel row.

### 13.3 AVX-512 backend

The AVX-512 backend:

- loads 64 bytes of packed BWT at a time,
- splits them into high and low nibbles,
- uses shuffle-based nibble lookup tables,
- reduces per-lane counts with `_mm512_sad_epu8`.

This backend is selected only if:

- AVX-512 support was compiled in,
- the host supports the relevant ISA and XCR0 state,
- `MAPPER_SIMD` resolves to `auto` or `avx512`.

### 13.4 Self-test coverage

`simd_selftest.cpp` checks:

- Occ correctness against a naive implementation,
- scalar/AVX-512 backend agreement when AVX-512 is available,
- multiple checkpoint strides and synthetic BWT patterns.

## 14. Benchmarking and Support Scripts

The project includes several support files:

- `bench.sh`
- `benchmark-memory.sh`
- `bench-mn5-subset.sh`
- `watch-progress.sh`
- `tail-progress.sh`
- `PROGRESS_MONITORING.md`
- `Makefile`

### 14.1 bench.sh

`bench.sh` is the main general-purpose benchmark harness.

It can:

- optionally build the index if missing,
- time mapper execution,
- validate the mapper output format,
- run a correctness comparison against `bwa mem`,
- record metrics and summaries.

Because the default index level is now `5`, the script now defaults `MAPPER_REF` to `REF`.

### 14.2 benchmark-memory.sh

This script sweeps levels `0` through `5` to compare:

- index size,
- index build time,
- mapper peak RSS,
- mapper runtime.

It is useful for development and memory-tradeoff exploration.

### 14.3 bench-mn5-subset.sh

This MN5-oriented helper runs:

- a subset benchmark,
- `MAPPER_SIMD=off`,
- `MAPPER_SIMD=avx512`,
- correctness checks against `bwa mem`,
- projected full-file runtime estimates.

It is intended for practical experimentation on MN5 where full single-thread runs can exceed the
available wall-clock budget.

### 14.4 Progress monitoring helpers

`watch-progress.sh` and `tail-progress.sh` are convenience scripts for following stderr progress
logs while a long benchmark is running. `PROGRESS_MONITORING.md` documents the intended workflow.

## 15. File-by-File Guide

This section maps the relevant repository files to their role.

### 15.1 Core headers

- `include/nucleotide.hpp`: DNA alphabet, normalization, reverse complement, and packed-storage helpers.
- `include/memory_stats.hpp`: current and peak RSS reporting.
- `include/bwt.hpp`: suffix-array/BWT public data structures.
- `include/fm_index.hpp`: FM-index configuration, owned/index-view structures, and serialization interface.
- `include/fm_search.hpp`: exact and inexact FM-search interface.
- `include/alignment.hpp`: banded alignment API and workspace.
- `include/fasta_reader.hpp`: sequential FASTA reader plus indexed streaming FASTA API.
- `include/fastq_reader.hpp`: FASTQ reader API.
- `include/simd_dispatch.hpp`: SIMD backend API and runtime-dispatch interface.

### 15.2 Core source files

- `src/indexer.cpp`: CLI for building FM-index files and selecting level `0..5`.
- `src/mapper.cpp`: CLI, threading, search orchestration, candidate verification, and output formatting.
- `src/bwt.cpp`: suffix-array construction and packed BWT generation.
- `src/fm_index.cpp`: FM-index building, serialization, mmap loading, Occ/LF/locate.
- `src/fm_search.cpp`: exact backward search and bounded inexact backtracking.
- `src/alignment.cpp`: banded edit-distance DP and simplified CIGAR traceback.
- `src/fasta_reader.cpp`: FASTA streaming and `.fai`-style indexed random access.
- `src/fastq_reader.cpp`: FASTQ parsing and validation.
- `src/simd_dispatch.cpp`: scalar backend and runtime backend selection.
- `src/simd_avx512.cpp`: AVX-512 packed-BWT counting backend.
- `src/simd_selftest.cpp`: regression and backend-agreement tests.

### 15.3 Build and benchmarking support

- `Makefile`: build modes, AVX-512 compilation, OpenMP, tests, and optional PGO/LTO.
- `bench.sh`: main benchmark harness with optional correctness comparison.
- `benchmark-memory.sh`: level sweep for memory/performance tradeoffs.
- `bench-mn5-subset.sh`: subset benchmark for MN5 with scalar vs AVX-512 comparison.
- `watch-progress.sh`: interactive progress dashboard for mapper stderr logs.
- `tail-progress.sh`: lighter real-time tail/grep helper.
- `PROGRESS_MONITORING.md`: notes about progress logging and debug builds.
- `README.md`: short project entry point.

### 15.4 Archived output directories

The directories:

- `benchmark-dataset-complet/`
- `benchmark-dataset-complet-max-level/`

are benchmark artifacts rather than source code. They are still useful because they capture real
measured summaries for the baseline and maximum-memory-optimization configurations.

## 16. Practical Recommendations

For normal use:

- if you care most about memory, use the default `L5` behavior,
- if you care most about speed, explicitly use `-L 0`,
- if you benchmark automatically after building an index, remember that `L5` requires `-R`.

For SIMD evaluation:

- compare `MAPPER_SIMD=off` and `MAPPER_SIMD=avx512` on the same subset,
- keep thread count fixed,
- do not assume AVX-512 helps equally across all levels,
- expect the effect to be larger when FM navigation dominates runtime.

For correctness:

- use the built-in `bench.sh` or `bench-mn5-subset.sh` comparison against `bwa mem`,
- treat `mapper_only`, `bwa_only`, and exact-position metrics as the key signals.

## 17. Summary

`mapper-memory` trades the large direct-address seed index for a chromosome-segmented FM-index.
That decision:

- dramatically reduces persistent index size,
- allows lower peak RSS during mapping,
- makes runtime much more sensitive to Occ/SA sparsity choices,
- and shifts the main performance bottleneck toward FM-index navigation.

The new default `L5` makes the project behave like a memory-contest mapper by default. The rest
of the codebase, scripts, and tests are now aligned with that assumption.

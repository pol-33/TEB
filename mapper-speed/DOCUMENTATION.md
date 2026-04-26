# mapper-speed Documentation

## 1. Overview

`mapper-speed` is the speed-oriented implementation of the final sequence-mapping project for the
course *Bioinformatics Algorithms and Techniques*. Its purpose is to align short sequencing reads
from a FASTQ file against a reference genome from a FASTA file, allowing a maximum edit
distance `k` with `0 <= k <= 3`, and to report at most the two best alignments for every read.

This implementation was designed with one main objective: maximize practical throughput while
remaining correct and staying within the rules of the project. In particular, the mapper is:

- single-threaded,
- optimized for short Illumina-like reads,
- optimized for repeated runs against the same reference,
- optimized for modern x86 CPUs with SIMD support,
- engineered to reject weak candidates as early and cheaply as possible.

The central idea is simple: most of the reference genome should never reach expensive dynamic
programming. The mapper therefore uses a seed-and-extend pipeline where cheap filters eliminate
nearly all bad candidates, and only a very small number of strong candidates reach the final
alignment stage.

## 2. Project Interface

The project supports the two-program interface described in the assignment.

### Index construction

```bash
./indexer -R genome.fa -I genome.compact.idx
```

### Mapping

```bash
./mapper -I genome.compact.idx -i reads.fastq -o output.sam -k 1
```

The `mapper` executable expects:

- `-I`: path to the index file,
- `-i`: path to the input FASTQ file,
- `-o`: path to the output SAM-like file,
- `-k`: maximum allowed edit distance, from `0` to `3`.

The `indexer` executable expects:

- `-R`: path to the reference FASTA file,
- `-I`: path where the index will be written.

Optional index construction modes:

```bash
./indexer -R genome.fa -I genome.dense.idx --mode dense
./indexer -R genome.fa -I genome.idx --mode auto
```

## 3. Recommended Mode for the Competition

The default indexing mode is `compact`, and this is the mode we recommend for the speed
competition.

This is not accidental. It was chosen after implementation and benchmarking because it provides
the best overall tradeoff between:

- runtime,
- memory usage,
- and mapping quality.

In other words:

- if the goal is to compete for **speed**, use `compact`,
- if the goal is to recover slightly more strong hits and memory is less of a concern, `dense`
  remains available.

So the recommended competition commands are:

```bash
./indexer -R genome.fa -I genome.compact.idx
./mapper -I genome.compact.idx -i reads.fastq -o output.sam -k 1
```

## 4. General Design Philosophy

The mapper is based on a seed-and-extend strategy.

Instead of scanning the whole genome with full edit-distance alignment, the pipeline is:

1. Build a fast index over the reference genome.
2. Choose a small set of informative seeds from the read.
3. Use seed hits to propose candidate mapping starts.
4. Reject weak candidates using packed exact-ish checks.
5. Apply a bounded verifier on promising equal-length candidates.
6. Run full banded alignment with traceback only for the best few candidates.
7. Report at most two best hits.

This architecture is what makes the mapper fast in practice. The final alignment stage is
important for correctness, but it is not used as the main search engine.

## 5. Reference Genome Representation

The reference-loading code is implemented in:

- `src/reference.cpp`
- `include/reference.hpp`

The reference is transformed into a compact internal representation:

- all bases are normalized,
- the complete genome is represented in a global coordinate space,
- standard bases are packed in **2 bits per base**,
- ambiguous bases (`N`) are stored in a separate bit-mask,
- chromosome boundaries and names are stored explicitly.

This representation is important for both speed and memory use.

Why it helps:

- packed storage reduces memory bandwidth,
- a separate `N` mask preserves correctness without forcing a slower mixed encoding,
- global coordinates simplify candidate generation and index lookup.

The index and the mapper therefore operate on one continuous genomic coordinate system, while
still being able to recover chromosome names and 1-based positions for output.

## 6. Index Structure

The index implementation lives mainly in:

- `src/indexer.cpp`
- `src/index.cpp`
- `include/index.hpp`

### 6.1 Seed model

The index is built on fixed-length seeds:

- seed length = `16`

Each valid 16-mer is encoded as a 32-bit integer key. That key identifies one logical entry in a
direct-address key space of size `4^16`.

This is conceptually very simple:

- every possible 16-mer has a unique integer key,
- the index maps that key to the list of reference positions where that seed occurs.

### 6.2 Why not a fully flat offset table?

In principle, a fully flat direct-address table is very fast, but it would be too large and wasteful.

To solve that, the index stores offsets in **compressed pages**:

- a page can be stored as a dense array when local variation is high,
- or as a sparse transition list when values are mostly constant.

This gives the best of both worlds:

- logically direct-address lookup,
- but much smaller memory footprint.

### 6.3 What is stored in the index?

The final on-disk index contains:

- index header,
- chromosome table,
- packed reference sequence,
- `N` mask,
- compressed offset pages,
- final seed-position array.

The index is memory-mapped (`mmap`) at runtime, so loading it is cheap and does not require
rebuilding large in-memory search structures from scratch.

## 7. Compact vs Dense Index

The project supports two index layouts.

### 7.1 Compact mode

Characteristics:

- default mode,
- stride = `4`,
- smaller index,
- lower RAM use,
- faster end-to-end mapping in our experiments,
- best practical mode for the benchmark competition.

Compact mode indexes fewer seed start positions in the reference. This slightly reduces candidate
coverage, but it improves cache behavior and reduces the amount of work during candidate
generation.

### 7.2 Dense mode

Characteristics:

- stride = `1`,
- larger index,
- higher RAM usage,
- slightly better sensitivity,
- slower than compact in practice.

Dense mode indexes every possible seed start, which gives richer coverage but also generates more
candidate work and more memory traffic.

### 7.3 Auto mode

`auto` estimates whether a dense index is likely to fit comfortably and chooses between dense and
compact accordingly. This mode is mainly useful as a convenience feature; for the competition, we
prefer making the choice explicit and using `compact`.

## 8. FASTQ Reading and Output Writing

The I/O code is implemented in:

- `src/buffered_io.cpp`

### 8.1 FASTQ reader

The FASTQ reader:

- reads sequentially using a large manual buffer,
- validates FASTQ structure,
- strips carriage returns when necessary,
- uses sequential-read hints on Linux,
- avoids per-line file I/O overhead.

This matters because the benchmark uses one million reads, so inefficient input handling would
quickly become visible.

### 8.2 Buffered output writer

The output writer:

- accumulates multiple records in memory,
- writes them in large chunks,
- avoids one system call per line.

This is especially important on HPC systems, where excessive small writes can become very
expensive.

## 9. Read Mapping Pipeline

The main search logic is implemented in:

- `src/search.cpp`
- `include/search.hpp`

For each input read, the mapper performs the following steps.

### 9.1 Read normalization and reverse-complement handling

The mapper considers both read orientations:

- original orientation,
- reverse complement.

Each orientation is processed independently, and the best hits from both orientations are merged
and ranked before final output.

### 9.2 Seed position selection

The mapper does not use every possible seed blindly. Instead, it chooses a small set of seed start
positions inside the read.

The policy depends on:

- read length,
- allowed edit distance `k`,
- and whether the index is dense or compact.

The goal is to balance two needs:

- cover the read well enough to find the correct region,
- but avoid too many seeds, because too many seed hits create too many candidates.

### 9.3 Seed extraction and frequency estimation

For the chosen read positions:

- rolling 16-mer keys are computed,
- invalid seeds are skipped,
- seed frequencies are obtained from the index,
- rare seeds are preferred.

This is a very important optimization. Frequent seeds tend to be uninformative and explode the
search space, especially in repetitive regions of the human genome.

### 9.4 Rare-seed-first candidate generation

Seed hits are converted into candidate starts by subtracting the read seed offset from the
reference hit position.

The candidate-generation stage includes several important details:

- chromosome-consistent start generation,
- limited `[-k, +k]` neighborhood expansion for indel tolerance,
- anchor sorting,
- duplicate anchor merging,
- candidate clustering,
- top-cluster selection with `partial_sort`.

The effect is that the mapper keeps only the most promising regions instead of trying to verify
everything.

## 10. Fast Early Acceptance Path

Before the full candidate machinery is used, the mapper tries a very cheap shortcut for easy reads.

If a few rare seeds produce a very small candidate set, the mapper can:

- compute packed mismatches immediately,
- accept exact or near-exact hits without entering the full search process.

This is particularly effective when:

- the read is unique,
- the true alignment is clean,
- and `k` is small.

On high-quality read data, this shortcut saves a meaningful amount of work.

## 11. Packed Mismatch Filtering

One of the most important speed optimizations in this mapper is that equal-length candidate
windows are first checked in packed form, before reconstructing the reference substring.

Implemented in:

- `src/index.cpp`
- `include/common.hpp`

The mapper compares:

- packed read bases,
- packed reference bases,
- and the reference `N` mask,

to count mismatches directly in compact form.

This stage is much cheaper than extracting a text substring and running full edit-distance
alignment immediately.

In practice, this filter removes many bad candidates before they become expensive.

## 12. Bounded Verifier

The verifier implementation is in:

- `src/verifier.cpp`
- `include/verifier.hpp`

The runtime supports dispatch between kernels depending on CPU features. However, an important
implementation decision was discovered during benchmarking:

- for the project setting `k <= 3`,
- the DP band width is only `2k+1 <= 7`,
- and forcing AVX512 into the verifier slowed the mapper down on MareNostrum 5.

Therefore, the current implementation intentionally prefers:

- `generic` on unsupported CPUs,
- `popcnt+bmi2` on modern x86 CPUs,

even when AVX512 features are available globally.

This is an example of a key engineering decision: not every theoretically more advanced SIMD path
is actually faster for the concrete workload.

## 13. Final Banded Alignment and CIGAR Recovery

The final alignment stage is implemented in:

- `src/alignment.cpp`

This stage runs only on a very small number of finalists and computes:

- edit distance,
- insertion/deletion/substitution path,
- final CIGAR string.

It uses:

- banded dynamic programming,
- traceback,
- CIGAR run-length compression.

This stage is the most expensive per-candidate operation, which is exactly why the rest of the
mapper is organized to keep the finalist set very small.

## 14. Simplified SAM-like Output

The output format follows the simplified project specification rather than the full SAM standard.

Each line contains:

- read name,
- chromosome name,
- 1-based mapping position,
- CIGAR,
- read sequence,
- quality string,
- optional `ALT:` field for the second-best hit.

If no valid hit exists, the mapper emits:

- `*` for chromosome,
- `0` for position,
- `*` for CIGAR.

The mapper reports at most two hits, always preferring the smallest edit distance.

## 15. CPU-Specific Optimizations

This mapper is single-threaded, but it is not scalar-only. It uses CPU-specific instructions where
they are useful and safe.

### 15.1 Feature detection

Detected CPU capabilities include:

- `popcnt`
- `bmi2`
- `avx`
- `avx2`
- `avx512f`
- `avx512dq`
- `avx512ifma`
- `avx512cd`
- `avx512bw`
- `avx512vl`
- `avx512vbmi`
- `avx512vbmi2`
- `avx512vnni`
- `avx512bitalg`
- `avx512vpopcntdq`
- `avx512fp16` when supported by the compiler

### 15.2 Where SIMD is really used

SIMD is used profitably in:

- bulk byte mismatch counting,
- packed mismatch counting,
- packed mismatch-word popcount,
- bulk compare helpers,
- other low-level counting primitives.

### 15.3 Why the verifier is still `popcnt+bmi2`

The verifier is the most subtle case.

For wider problems, AVX512 can be useful. But for this project:

- `k` is tiny,
- the active band is very narrow,
- and the cost of AVX512 setup and possible CPU downclock can exceed the benefit.

So the optimized decision is not “use AVX512 everywhere”, but “use it only where it helps”.

This is why the benchmark may show:

- `effective = avx512+...`
- but `verifier = popcnt+bmi2`

and that is the intended behavior.

## 16. Memory and Cache Optimizations

A large part of the mapper’s speed comes from careful memory engineering.

Important optimizations include:

- packed 2-bit reference storage,
- separate `N` bit-mask,
- page-compressed offset lookup,
- memory-mapped index loading,
- chromosome lookup bins,
- scratch-buffer reuse,
- in-place anchor merging,
- `partial_sort` instead of full sort when only top results matter,
- buffered output,
- software prefetching of:
  - occurrence-list entries,
  - future reference windows,
  - packed-reference and `N`-mask regions.

These are not “nice extras”; they are essential to practical performance on genome-scale data.

## 17. Build-Level Optimizations

The build system in `Makefile` was tuned for optimized release binaries.

Current optimized configuration includes:

- `-O3`
- `-DNDEBUG`
- `-fomit-frame-pointer`
- `-fstrict-aliasing`
- `-march=native`
- `-mtune=native`
- `-flto` on Linux
- `-fno-semantic-interposition` on Linux

Additionally, the project supports optional PGO:

```bash
make pgo-generate
# run a representative workload
make pgo-use
```

PGO is useful because it allows the compiler to specialize the binary for the code paths that are
actually frequent in the benchmark workload.

## 18. Correct Execution on MareNostrum 5

One of the practical lessons learned during benchmarking is that running directly from the shared
filesystem can make the mapper appear slow or even apparently “stuck”, even when the algorithm
is correct.

The correct way to benchmark on MN5 is:

- copy the reference,
- copy the reads,
- copy the mapper indexes,
- copy the prebuilt BWA index,
- execute from `TMPDIR`.

Recommended benchmark pattern:

```bash
mkdir -p "$TMPDIR/mapper-bench" && \
cp ~/TEB/mapper-speed/genome.dense.idx "$TMPDIR/mapper-bench/" && \
cp ~/TEB/mapper-speed/genome.compact.idx "$TMPDIR/mapper-bench/" && \
cp ~/TEB/data/reads_1M.fastq "$TMPDIR/mapper-bench/" && \
cp ~/TEB/data/GRCh38.fna "$TMPDIR/mapper-bench/" && \
mkdir -p "$TMPDIR/mapper-bench/bwa-index" && \
cp ~/TEB/mapper-speed/bwa-index/reference.bwa.* "$TMPDIR/mapper-bench/bwa-index/" && \
cd ~/TEB/mapper-speed && \
DENSE_INDEX="$TMPDIR/mapper-bench/genome.dense.idx" \
COMPACT_INDEX="$TMPDIR/mapper-bench/genome.compact.idx" \
REF="$TMPDIR/mapper-bench/GRCh38.fna" \
READS="$TMPDIR/mapper-bench/reads_1M.fastq" \
OUT_DIR="$TMPDIR/mapper-bench/out" \
BWA_PREFIX="$TMPDIR/mapper-bench/bwa-index/reference.bwa" \
BUILD=0 \
BWA_THREADS=1 \
./bench-mn5.sh
```

This avoids shared-filesystem stalls and produces much more stable performance measurements.

## 19. Benchmark Results Observed on MN5

One representative benchmark on MareNostrum 5 gave the following results for `1M` reads and
`k = 1`:

### Compact native

- runtime: `91.52 s`
- peak RSS: `7,274,393,600 bytes` (`~6.77 GiB`)
- mapped reads: `804,366`
- unmapped reads: `195,634`
- reads with alternative hit: `257,916`

Correctness versus BWA:

- shared mapped reads: `803,620`
- same chromosome: `766,278`
- within 10 bp: `566,821`
- exact position: `566,819`
- exact CIGAR: `765,912`
- exact primary agreement: `566,818`

### Dense native

- runtime: `122.72 s`
- peak RSS: `19,949,342,720 bytes` (`~18.58 GiB`)
- mapped reads: `806,658`
- unmapped reads: `193,342`
- reads with alternative hit: `233,880`

Correctness versus BWA:

- shared mapped reads: `805,885`
- same chromosome: `775,218`
- within 10 bp: `596,623`
- exact position: `596,623`
- exact CIGAR: `774,939`
- exact primary agreement: `596,622`

### BWA baseline

- runtime: `330.65 s`
- peak RSS: `6,129,537,024 bytes` (`~5.71 GiB`)

### Interpretation

These measurements show:

- `compact` is the fastest mode and the best competition choice,
- `dense` is slower and much heavier in memory, but slightly stronger in hit recovery,
- both mapper modes are much faster than `bwa mem` in this benchmark,
- `compact-native` is roughly `3.6x` faster than BWA,
- `dense-native` is roughly `2.7x` faster than BWA.

This is exactly why `compact` is the default and the recommended mode for the speed competition.

## 20. Why the Mapper Is Fast

The performance comes from combining many compatible decisions:

- build an index once and reuse it,
- use compact mode by default,
- represent the reference in packed form,
- use direct-address logical seed lookup,
- choose rare informative seeds first,
- prune candidate regions aggressively,
- exploit a fast exact-ish shortcut,
- compare packed reference windows before extracting text,
- use a bounded verifier before full alignment,
- align only the best few candidates,
- exploit useful CPU instructions without forcing harmful ones,
- optimize memory traffic and cache behavior,
- compile with aggressive release flags.

There is no single trick that explains the speed. The mapper is fast because all stages were
designed to avoid unnecessary work.

## 21. Final Recommendation

For the course project and especially for the speed competition, the recommended workflow is:

1. Build the compact index:

```bash
./indexer -R genome.fa -I genome.compact.idx
```

2. Run the mapper:

```bash
./mapper -I genome.compact.idx -i reads.fastq -o output.sam -k 1
```

3. On MN5, run from `TMPDIR`.

4. If extra tuning is desired, build with PGO after obtaining representative profiling runs.

In summary:

- use `compact` for speed,
- use `dense` only when slightly higher candidate coverage is worth the extra cost,
- keep the implementation single-threaded,
- and rely on the current filtering and packed-data pipeline to achieve high practical throughput.

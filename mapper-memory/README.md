# mapper-memory

`mapper-memory` is the memory-oriented mapper implementation for the TEB sequence-mapper project.

## What It Builds

- `indexer`: builds a serialized FM-index from a reference FASTA
- `mapper`: maps FASTQ reads using the serialized FM-index via `mmap`

## Build

```bash
make
```

## Commands

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.idx
./mapper  -I /path/to/genome.idx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
./mapper  -R /path/to/genome.fa  -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

`-I` mode is the intended competition workflow. `-R` mode is a correctness-first fallback that scans the genome directly and is not designed for the contest-sized input.

## Benchmarking

Use the bundled benchmark harness to run a timed mapper job, validate the simplified SAM-like output, and optionally compare a smaller subset against `bwa mem`:

```bash
./bench.sh
```

Useful overrides:

```bash
BUILD_INDEX_IF_MISSING=1 BENCH_READS=10000 CORRECTNESS_READS=500 ./bench.sh
BUILD_INDEX_IF_MISSING=1 BENCH_READS=0 ./bench.sh
BUILD_BWA_INDEX_IF_MISSING=1 ./bench.sh
```

By default the script avoids surprise multi-hour setup work: it will reuse an existing FM-index, and it will only build the FM-index or the `bwa` index automatically when the corresponding `*_IF_MISSING=1` flag is set.

## Implementation Notes

- FASTA input is flattened into a contiguous genome plus chromosome metadata.
- The suffix array is built with SA-IS, then converted into a packed BWT.
- The serialized index stores the primary row, sampled `Occ`, sampled `SA`, chromosome metadata, and a packed copy of the reference sequence so the mapper can extract local windows for CIGAR generation without re-reading the FASTA.
- The mapper streams FASTQ records, performs bounded FM-index backtracking for `k <= 3`, verifies candidates with banded edit-distance DP, and writes up to two alignments per read in the simplified SAM-like contest format.

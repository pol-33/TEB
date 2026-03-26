# mapper-memory

`mapper-memory` is the memory-oriented mapper implementation for the TEB sequence-mapper project.

## What It Builds

- `indexer`: builds a serialized low-memory seed index from a reference FASTA
- `mapper`: maps FASTQ reads against the serialized seed index via `mmap`

## Build

```bash
make
```

## Commands

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.seed.idx
./mapper  -I /path/to/genome.seed.idx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

`-I` mode is the intended competition workflow. The contest-scale mapper currently requires a prebuilt index.

## Benchmarking

Use the bundled benchmark harness to run a timed mapper job, validate the simplified SAM-like output, and optionally compare a smaller subset against `bwa mem`:

```bash
./bench.sh
```

Useful overrides:

```bash
BUILD_INDEX_IF_MISSING=1 BENCH_READS=10000 CORRECTNESS_READS=500 ./bench.sh
INDEX=./bench-results/genome.seed.idx BENCH_READS=500 RUN_BWA_CHECK=0 ./bench.sh
BUILD_BWA_INDEX_IF_MISSING=1 ./bench.sh
```

By default the script avoids surprise multi-hour setup work: it will reuse an existing seed index, and it will only build the seed index or the `bwa` index automatically when the corresponding `*_IF_MISSING=1` flag is set.

## Implementation Notes

- FASTA input is flattened into a contiguous genome plus chromosome metadata.
- The on-disk index stores a packed 2-bit reference plus a direct-address seed table for 13-mers.
- Index construction uses bounded partitioned writes so GRCh38 can be indexed without the old suffix-array memory failure.
- The mapper streams FASTQ records, searches both forward and reverse-complement orientations, verifies candidates with banded edit-distance DP, and writes up to two alignments per read in the simplified SAM-like contest format.

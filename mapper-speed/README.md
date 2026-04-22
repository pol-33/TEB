# mapper-speed

`mapper-speed` is the speed-focused implementation for the TEB mapper project.

## Design

- Whole-reference `16-mer` seed index with direct-address logical lookup.
- Global coordinate space across all chromosomes.
- Candidate voting from rare seeds first.
- Two-stage verification:
  - bounded Myers edit-distance filter,
  - banded DP with CIGAR recovery for finalists.
- Single-thread mapping, tuned for short Illumina-like reads and `k <= 3`.

## Build

```bash
make
```

Portable build:

```bash
make portable
```

## Commands

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.idx
./mapper  -I /path/to/genome.idx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

## Smoke Test

```bash
make test
```

## Benchmark

```bash
BENCH_READS=1000 ./bench.sh
```

MareNostrum 5 four-way benchmark:

```bash
sbatch mapper-speed/bench-mn5.sh
```

Useful overrides:

```bash
DENSE_INDEX=$TMPDIR/genome.idx \
COMPACT_INDEX=$TMPDIR/genome.compact.idx \
READS=$TMPDIR/reads_1M.fastq \
OUT_DIR=$PWD/mapper-speed/bench-mn5-results \
sbatch mapper-speed/bench-mn5.sh
```

To disable AVX512 for a manual run, cap the runtime dispatch:

```bash
MAPPER_SPEED_MAX_SIMD=avx2 ./mapper -I genome.idx -i reads.fastq -o out.sam -k 1
```

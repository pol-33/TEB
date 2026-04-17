# mapper-memory

`mapper-memory` is the memory-contest mapper for the TEB sequence-mapper project.

## What It Builds

- `indexer`: builds a serialized FM-index from a reference FASTA
- `mapper`: maps FASTQ reads against the memory-mapped FM-index

## Build

```bash
make
```

Optional SIMD modes:

```bash
make SIMD=auto    # default: portable build, use AVX-512 only when supported
make SIMD=off     # scalar-only build
make SIMD=avx512  # require x86_64 compiler support for AVX-512
```

Run the tiny smoke test:

```bash
make test
```

Run the benchmark harness:

```bash
./bench.sh
```

The harness auto-detects whether `/usr/bin/time -l` is available. On restricted macOS environments it falls back to `/usr/bin/time -p` and uses the mapper/indexer peak-RSS report from stderr.

## Commands

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.fmidx
./mapper  -I /path/to/genome.fmidx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

Optional runtime override for the packed-BWT occurrence counter:

```bash
MAPPER_SIMD=auto ./mapper ...
MAPPER_SIMD=off ./mapper ...
MAPPER_SIMD=avx512 ./mapper ...
```

The contest-scale mapper currently requires a prebuilt index. Reads are streamed one record at a time, and the output format is:

```text
read_name chrom pos_1based cigar seq qual [ALT:chrom,pos,cigar]
```

Unmapped reads are emitted as:

```text
read_name * 0 * seq qual
```

## Implementation Notes

- The indexer builds one native FM-index segment per chromosome, which keeps suffix-array construction bounded by chromosome size instead of the full reference size.
- Each chromosome stores its own packed BWT, Occ checkpoints, sampled SA, and packed reference sequence inside one mmap-friendly index file.
- The mapper uses exact backward search first, then bounded inexact FM backtracking for `k <= 3`.
- Candidate locations are verified with the reusable banded edit-distance aligner and written as the best two simplified SAM-like hits.
- The mapper searches both forward and reverse-complement orientations.

## Memory Analysis

Projected GRCh38 FM-index footprint:

| Component | Size (GRCh38) | RSS contribution (mmap) |
| --- | ---: | ---: |
| Packed BWT (2b/base) | ~0.775 GB | ~0.2-0.5 GB |
| Occ checkpoints (256 stride) | ~0.193 GB | ~0.1-0.2 GB |
| Sampled SA (1/32) | ~0.387 GB | ~0.1-0.2 GB |
| Packed genome | ~0.775 GB | ~0.05-0.2 GB |
| Total index on disk | ~2.13 GB | ~0.5-1.1 GB peak RSS (estimated) |

Comparison with the old seed-index approach:

| Component | Size |
| --- | ---: |
| Positions array | ~12.4 GB |
| Bucket offsets | ~0.27 GB |
| Total | ~13.4 GB |

The FM-index wins the memory contest because the mapper touches only the BWT, Occ checkpoints, sampled SA, and small reference windows through `mmap`, instead of paging in explicit position lists for nearly every seed in the genome.

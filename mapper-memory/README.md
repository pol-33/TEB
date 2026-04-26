# mapper-memory

`mapper-memory` is the memory-focused implementation for the TEB mapper project.

## Authors

- Pol Plana Torrents
- Joan Teruel
- Mariam Delgado

## Documentation

For the full project description, file-by-file design notes, indexing-level tradeoffs,
SIMD dispatch details, and benchmarking workflow, see [DOCUMENTATION.md](DOCUMENTATION.md).

## Build

```bash
make
```

Portable build:

```bash
make portable
```

## Quick Start

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.fmidx
./mapper  -I /path/to/genome.fmidx -i /path/to/reads.fastq -o /path/to/output.sam -k 1 -R /path/to/genome.fa
```

`indexer` now defaults to optimization level `5`, which minimizes index size and RAM usage.
Because level `5` does not embed the packed reference genome, the default `mapper` command
should also receive `-R /path/to/genome.fa`.

If you want the old largest/fastest baseline index with the genome embedded inside the FM-index:

```bash
./indexer -R /path/to/genome.fa -I /path/to/genome.fmidx -L 0
./mapper  -I /path/to/genome.fmidx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

## Smoke Test

```bash
make test
```

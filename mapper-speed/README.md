# mapper-speed

`mapper-speed` is the speed-focused implementation for the TEB mapper project.

## Authors

- Pol Plana Torrents
- Joan Teruel
- Mariam Delgado

## Documentation

For the full project description, design rationale, indexing modes, benchmark workflow,
profiling workflow, and MN5 notes, see [DOCUMENTATION.md](DOCUMENTATION.md).

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
./indexer -R /path/to/genome.fa -I /path/to/genome.compact.idx
./mapper  -I /path/to/genome.compact.idx -i /path/to/reads.fastq -o /path/to/output.sam -k 1
```

`indexer` uses the compact index layout by default. We recommend using the compact index layout for all applications, as it is faster to load and has a smaller memory footprint.

## Smoke Test

```bash
make test
```

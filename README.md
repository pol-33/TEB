# TEB Sequence Mapper Project

This repository now contains two mapper variants for the TEB competitive bioinformatics project:

- `mapper-memory/`: a low-memory seed-index mapper tuned for lower peak RSS during mapping
- `mapper-speed/`: a future minimizer/hash mapper reserved for the speed-focused track

The original FASTA/FASTQ parser coursework has been preserved under `legacy/`.

## Build

Build the memory-oriented implementation:

```bash
cd mapper-memory
make
```

Build the speed-oriented placeholder:

```bash
cd mapper-speed
make
```

## Data

Recommended inputs:

- GRCh38 reference FASTA as `genome.fa`
- First 1 million reads from ERR194147 as `reads_1M.fastq`

Suggested download commands:

```bash
wget -O genome.fa https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip genome.fa
prefetch ERR194147
fasterq-dump ERR194147 --split-files
head -n 4000000 ERR194147_1.fastq > reads_1M.fastq
```

See [data/README.md](/Users/polplana/Documents/Universitat/TEB/teb/data/README.md) for local placement guidance.

## CLI

Build the memory mapper index:

```bash
cd mapper-memory
./indexer -R /path/to/genome.fa -I /path/to/genome.seed.idx
```

Map reads with the serialized index:

```bash
cd mapper-memory
./mapper -I /path/to/genome.seed.idx -i /path/to/reads_1M.fastq -o /path/to/output.sam -k 1
```

Output format:

```text
read_name chrom pos_1based cigar seq qual [ALT:chrom,pos,cigar]
```

Unmapped reads are emitted as:

```text
read_name * 0 * seq qual
```

## Expected Performance

Approximate target numbers for the memory-oriented mapper on contest-sized data:

- Index build: a large one-time job on GRCh38
- Mapping: depends strongly on `k`, read repetitiveness, and file-cache warmth

These numbers are project goals rather than guarantees for every machine. The current implementation prioritizes the required architecture and correctness path first.

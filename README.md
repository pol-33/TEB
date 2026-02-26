# TEB — Bioinformatics Parser

## Build

```bash
make        # builds the `teb` executable
make clean  # removes object files and the executable
make re     # clean + build
```

The exercises can be compiled independently:
```bash
make exercises/ex2_1   # naive exact search
make exercises/ex5_2   # full Boyer-Moore
```

---

## Usage

```
./teb --input <file> --format <fasta|fastq> [--output <file>] [--qmin <int>]
```

| Argument | Required | Description |
|---|---|---|
| `--input <file>` | ✅ | Input FASTA or FASTQ file |
| `--format <fasta\|fastq>` | ✅ | Input file format |
| `--output <file>` | ❌ | Write parsed sequences to this file. Omit for stats-only mode |
| `--qmin <int>` | ❌ | FASTQ only — trim bases from the right with quality score below this threshold |

---

## Examples

**FASTA — stats only (no output file):**
```bash
./teb --input datasets/chr1.fna --format fasta
```

**FASTA — parse and write output:**
```bash
./teb --input datasets/chr1.fna --format fasta --output output/chr1-out.fna
```

**FASTQ — stats only:**
```bash
./teb --input datasets/SRR22320000_1.fastq --format fastq
```

**FASTQ — parse with quality trimming and write output:**
```bash
./teb --input datasets/SRR22320000_1.fastq --format fastq --output output/SRR1-out.fastq --qmin 20
```

---

## Project structure

```
teb/
├── src/               # Source files (.cpp)
├── include/           # Header files (.hpp)
├── exercises/         # Standalone exercise programs
├── datasets/          # Input genomic data
├── docs/              # Lab assignment PDFs
├── output/            # Generated output files
└── Makefile
```

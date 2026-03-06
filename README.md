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
./teb.exe -i <file> -f <fasta|fastq> [-o <file>] [-qmin <int>]
```

| Argument | Required | Description |
|---|---|---|
| `-i <file>` | ✅ | Input FASTA or FASTQ file |
| `-f <fasta\|fastq>` | ✅ | Input file format |
| `-o <file>` | ❌ | Write parsed sequences to this file. Omit for stats-only mode |
| `-qmin <int>` | ❌ | FASTQ only — trim bases from the right with quality score below this threshold |
| `-k <int>` | ❌ | FASTA only — compute k-mer frequencies for this k (e.g. 3 for trinucleotides) |

---

## Examples

**FASTA — stats only (no output file):**
```bash
./teb.exe -i datasets/chr1.fna -f fasta
```

**FASTA — parse and write output:**
```bash
./teb.exe -i datasets/chr1.fna -f fasta -o output/chr1-out.fna
```

**FASTQ — stats only:**
```bash
./teb.exe -i datasets/SRR22320000_1.fastq -f fastq
```

**FASTQ — parse with quality trimming and write output:**
```bash
./teb.exe -i datasets/SRR22320000_1.fastq -f fastq -o output/SRR1-out.fastq -qmin 20
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

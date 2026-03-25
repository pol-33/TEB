# Data Layout

The competition workflow expects two main files:

- `genome.fa`: the GRCh38 reference FASTA
- `reads_1M.fastq`: the first 1 million 100 bp reads derived from ERR194147

You can place them wherever you want and point `indexer` / `mapper` at their full paths. If you want a simple local convention, place them under `data/`:

```text
data/
├── genome.fa
└── reads_1M.fastq
```

This repository also keeps the pre-existing `datasets/` directory at the top level for legacy parser inputs from the original coursework project.

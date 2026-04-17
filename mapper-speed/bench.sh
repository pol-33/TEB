#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$ROOT_DIR/.." && pwd)"

REF="${REF:-$REPO_DIR/data/genome.fa}"
READS="${READS:-$REPO_DIR/data/reads_1M.fastq}"
INDEX="${INDEX:-$ROOT_DIR/bench-results/genome.speed.idx}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/bench-results}"
K="${K:-1}"
BENCH_READS="${BENCH_READS:-1000}"

mkdir -p "$OUT_DIR"

SUBSET="$OUT_DIR/reads_subset.fastq"
if [[ "$BENCH_READS" -gt 0 ]]; then
  awk -v limit="$((BENCH_READS * 4))" 'NR <= limit { print }' "$READS" > "$SUBSET"
else
  cp "$READS" "$SUBSET"
fi

make -C "$ROOT_DIR" -j4

if [[ ! -f "$INDEX" ]]; then
  echo "[bench] building index"
  /usr/bin/time -p -o "$OUT_DIR/index.time.log" "$ROOT_DIR/indexer" -R "$REF" -I "$INDEX" \
    > "$OUT_DIR/index.stdout.log" 2> "$OUT_DIR/index.stderr.log"
fi

echo "[bench] running mapper"
/usr/bin/time -p -o "$OUT_DIR/mapper.time.log" "$ROOT_DIR/mapper" -I "$INDEX" -i "$SUBSET" -o "$OUT_DIR/out.sam" -k "$K" \
  > "$OUT_DIR/mapper.stdout.log" 2> "$OUT_DIR/mapper.stderr.log"

REAL_SECONDS="$(awk '/^real / {print $2}' "$OUT_DIR/mapper.time.log")"
READ_COUNT="$(awk 'END { printf "%d\n", NR / 4 }' "$SUBSET")"
RPS="$(awk -v reads="$READ_COUNT" -v real="$REAL_SECONDS" 'BEGIN { if (real > 0) printf "%.2f", reads / real; else print "0.00" }')"
MAPPED="$(awk '$2 != "*" { count++ } END { print count + 0 }' "$OUT_DIR/out.sam")"
UNMAPPED="$(awk '$2 == "*" { count++ } END { print count + 0 }' "$OUT_DIR/out.sam")"

{
  echo "Benchmark Summary"
  echo "================="
  echo "reference: $REF"
  echo "reads: $READS"
  echo "bench_reads: $READ_COUNT"
  echo "index: $INDEX"
  echo "k: $K"
  echo
  echo "Mapping Metrics"
  echo "---------------"
  echo "real_seconds: $REAL_SECONDS"
  echo "reads_per_second: $RPS"
  echo "mapped: $MAPPED"
  echo "unmapped: $UNMAPPED"
} | tee "$OUT_DIR/summary.txt"

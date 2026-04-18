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
COMPARE_BWA="${COMPARE_BWA:-0}"
RUNS="${RUNS:-1}"
BWA_THREADS="${BWA_THREADS:-1}"
BWA_PREFIX="${BWA_PREFIX:-$OUT_DIR/reference.bwa}"

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

run_avg() {
  local runs="$1"
  shift
  local sum=0
  for i in $(seq 1 "$runs"); do
    local real
    real="$(
      TIMEFORMAT=%R
      { time "$@" >/dev/null 2>/dev/null; } 2>&1
    )"
    sum="$(awk -v a="$sum" -v b="$real" 'BEGIN { printf "%.6f", a + b }')"
  done
  awk -v s="$sum" -v n="$runs" 'BEGIN { if (n > 0) printf "%.6f", s / n; else print "0.000000" }'
}

echo "[bench] running mapper"
MAPPER_AVG_REAL="$(run_avg "$RUNS" "$ROOT_DIR/mapper" -I "$INDEX" -i "$SUBSET" -o "$OUT_DIR/out.sam" -k "$K")"
"$ROOT_DIR/mapper" -I "$INDEX" -i "$SUBSET" -o "$OUT_DIR/out.sam" -k "$K" \
  > "$OUT_DIR/mapper.stdout.log" 2> "$OUT_DIR/mapper.stderr.log"

READ_COUNT="$(awk 'END { printf "%d\n", NR / 4 }' "$SUBSET")"
RPS="$(awk -v reads="$READ_COUNT" -v real="$MAPPER_AVG_REAL" 'BEGIN { if (real > 0) printf "%.2f", reads / real; else print "0.00" }')"
MAPPED="$(awk '$2 != "*" { count++ } END { print count + 0 }' "$OUT_DIR/out.sam")"
UNMAPPED="$(awk '$2 == "*" { count++ } END { print count + 0 }' "$OUT_DIR/out.sam")"
MAPPER_PEAK_RSS="$(awk '/peak RSS:/ {print $(NF-1)}' "$OUT_DIR/mapper.stderr.log" | tail -n 1)"

if [[ "$COMPARE_BWA" == "1" ]]; then
  if [[ ! -f "${BWA_PREFIX}.bwt" ]]; then
    echo "[bench] building bwa index"
    /usr/bin/time -p -o "$OUT_DIR/bwa.index.time.log" bwa index -p "$BWA_PREFIX" "$REF" \
      > "$OUT_DIR/bwa.index.stdout.log" 2> "$OUT_DIR/bwa.index.stderr.log"
  fi
  echo "[bench] running bwa mem"
  BWA_AVG_REAL="$(run_avg "$RUNS" bwa mem -t "$BWA_THREADS" "$BWA_PREFIX" "$SUBSET")"
  bwa mem -t "$BWA_THREADS" "$BWA_PREFIX" "$SUBSET" \
    > "$OUT_DIR/bwa.sam" 2> "$OUT_DIR/bwa.stderr.log"
  BWA_RPS="$(awk -v reads="$READ_COUNT" -v real="$BWA_AVG_REAL" 'BEGIN { if (real > 0) printf "%.2f", reads / real; else print "0.00" }')"
  BWA_MAPPED="$(awk 'BEGIN { FS="\t" } /^@/ { next } $3 != "*" { count++ } END { print count + 0 }' "$OUT_DIR/bwa.sam")"
  BWA_UNMAPPED="$(awk 'BEGIN { FS="\t" } /^@/ { next } $3 == "*" { count++ } END { print count + 0 }' "$OUT_DIR/bwa.sam")"
fi

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
  echo "mapper_real_seconds_avg: $MAPPER_AVG_REAL"
  echo "mapper_reads_per_second: $RPS"
  echo "mapper_peak_rss_mib: ${MAPPER_PEAK_RSS:-unknown}"
  echo "mapper_mapped: $MAPPED"
  echo "mapper_unmapped: $UNMAPPED"
  if [[ "$COMPARE_BWA" == "1" ]]; then
    echo
    echo "bwa mem Metrics"
    echo "---------------"
    echo "bwa_real_seconds_avg: $BWA_AVG_REAL"
    echo "bwa_reads_per_second: $BWA_RPS"
    echo "bwa_mapped: $BWA_MAPPED"
    echo "bwa_unmapped: $BWA_UNMAPPED"
    echo "runs: $RUNS"
    echo "threads: mapper=1 bwa=$BWA_THREADS"
  fi
} | tee "$OUT_DIR/summary.txt"

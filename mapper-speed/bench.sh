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
RUNS="${RUNS:-1}"
COMPARE_BWA="${COMPARE_BWA:-1}"
BWA_THREADS="${BWA_THREADS:-1}"
BWA_BIN="${BWA_BIN:-bwa}"
BWA_PREFIX="${BWA_PREFIX:-$OUT_DIR/reference.bwa}"
TIME_STYLE="${TIME_STYLE:-auto}"

mkdir -p "$OUT_DIR"

log() {
  printf '[bench] %s\n' "$*"
}

die() {
  printf '[bench][error] %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: ./bench.sh

Environment overrides:
  REF=...             Reference FASTA
  READS=...           Input FASTQ
  INDEX=...           mapper-speed index
  OUT_DIR=...         Output directory for logs and summaries
  K=1                 Maximum edit distance
  BENCH_READS=1000    Reads to benchmark; set to 0 for all reads
  RUNS=1              Number of timing runs per tool
  COMPARE_BWA=1       Set to 0 to benchmark only mapper-speed
  BWA_THREADS=1       Threads for bwa mem
  BWA_BIN=bwa         Path to bwa executable
  BWA_PREFIX=...      bwa index prefix path
  TIME_STYLE=auto     auto | full | portable
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ -f "$REF" ]] || die "reference FASTA not found: $REF"
[[ -f "$READS" ]] || die "reads FASTQ not found: $READS"
[[ -f "$INDEX" ]] || die "mapper index not found: $INDEX"
[[ "$K" =~ ^[0-9]+$ ]] || die "K must be an integer"
[[ "$RUNS" =~ ^[0-9]+$ ]] || die "RUNS must be an integer"
[[ "$BENCH_READS" =~ ^[0-9]+$ ]] || die "BENCH_READS must be an integer"

count_reads() {
  awk 'END { printf "%d\n", NR / 4 }' "$1"
}

extract_subset() {
  local input_path="$1"
  local read_count="$2"
  local output_path="$3"
  if [[ "$read_count" -le 0 ]]; then
    cp "$input_path" "$output_path"
  else
    awk -v limit="$((read_count * 4))" 'NR <= limit { print }' "$input_path" > "$output_path"
  fi
}

detect_time_style() {
  case "$TIME_STYLE" in
    full|portable)
      printf '%s\n' "$TIME_STYLE"
      return
      ;;
    auto)
      ;;
    *)
      die "TIME_STYLE must be one of: auto, full, portable"
      ;;
  esac

  local probe
  probe="$(mktemp)"
  if /usr/bin/time -l -o "$probe" true >/dev/null 2>&1; then
    rm -f "$probe"
    printf 'full\n'
  else
    rm -f "$probe"
    printf 'portable\n'
  fi
}

parse_time_report() {
  local time_path="$1"
  awk '
    /real/ && /user/ && /sys/            { real=$1; user=$3; sys=$5 }
    /^real[[:space:]]+/                  { real=$2 }
    /^user[[:space:]]+/                  { user=$2 }
    /^sys[[:space:]]+/                   { sys=$2 }
    /maximum resident set size/          { rss=$1 }
    /page faults/                        { pf=$1 }
    END {
      if (real == "") real = 0;
      if (user == "") user = 0;
      if (sys == "") sys = 0;
      if (rss == "") rss = 0;
      if (pf == "") pf = 0;
      printf "real=%s\nuser=%s\nsys=%s\nrss=%s\npf=%s\n", real, user, sys, rss, pf;
    }
  ' "$time_path"
}

parse_mapper_peak_rss_bytes() {
  local stderr_path="$1"
  awk '
    /peak RSS:/ {
      line = $0;
      sub(/^.*peak RSS: /, "", line);
      sub(/ MiB.*$/, "", line);
      rss = int((line * 1024 * 1024) + 0.5);
    }
    END {
      if (rss == "") rss = 0;
      print rss;
    }
  ' "$stderr_path"
}

run_once() {
  local label="$1"
  local stdout_path="$2"
  local stderr_path="$3"
  local time_path="$4"
  shift 4

  if [[ "$RESOLVED_TIME_STYLE" == "full" ]]; then
    /usr/bin/time -l -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
  else
    /usr/bin/time -p -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
  fi
}

run_avg() {
  local runs="$1"
  shift
  local sum=0
  local i
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

count_mapper_stats() {
  local sam_path="$1"
  local report_path="$2"
  awk '
    BEGIN { FS = "[ \t]+" }
    {
      reads++;
      if ($2 == "*") {
        unmapped++;
      } else {
        mapped++;
      }
      if (NF >= 7 && $7 ~ /^ALT:/) {
        with_alt++;
      }
    }
    END {
      printf "reads=%d\nmapped=%d\nunmapped=%d\nwith_alt=%d\n",
             reads + 0, mapped + 0, unmapped + 0, with_alt + 0;
    }
  ' "$sam_path" > "$report_path"
}

count_bwa_stats() {
  local sam_path="$1"
  local max_nm="$2"
  local report_path="$3"
  awk -v max_nm="$max_nm" '
    function bit(flag, mask) { return int(flag / mask) % 2; }
    function valid_cigar(cigar) { return cigar ~ /^([0-9]+[MID])+$/; }
    /^@/ { next }
    {
      flag = $2 + 0;
      name = $1;
      is_secondary = bit(flag, 256) || bit(flag, 2048);
      is_unmapped = ($3 == "*" || bit(flag, 4));
      nm = -1;
      for (i = 12; i <= NF; ++i) {
        if ($i ~ /^NM:i:/) {
          nm = substr($i, 6) + 0;
        }
      }
      is_valid = (!is_unmapped && nm >= 0 && nm <= max_nm && valid_cigar($6));
      if (!is_secondary) {
        reads++;
        if (is_valid) mapped++;
        else unmapped++;
      } else if (is_valid) {
        alt[name] = 1;
      }
    }
    END {
      for (name in alt) with_alt++;
      printf "reads=%d\nmapped=%d\nunmapped=%d\nwith_alt=%d\n",
             reads + 0, mapped + 0, unmapped + 0, with_alt + 0;
    }
  ' "$sam_path" > "$report_path"
}

convert_mapper_output() {
  local input_path="$1"
  local output_path="$2"
  awk 'BEGIN { FS = "[ \t]+" } { print $1 "\t" $2 "\t" $3 "\t" $4 }' "$input_path" > "$output_path"
}

convert_bwa_output() {
  local input_path="$1"
  local output_path="$2"
  local max_nm="$3"
  awk -v max_nm="$max_nm" '
    function bit(flag, mask) { return int(flag / mask) % 2; }
    function valid_cigar(cigar) { return cigar ~ /^([0-9]+[MID])+$/; }
    /^@/ { next }
    {
      flag = $2 + 0;
      if (bit(flag, 256) || bit(flag, 2048)) {
        next;
      }
      if ($3 == "*" || bit(flag, 4)) {
        print $1 "\t*\t0\t*";
        next;
      }
      nm = -1;
      for (i = 12; i <= NF; ++i) {
        if ($i ~ /^NM:i:/) {
          nm = substr($i, 6) + 0;
        }
      }
      if (nm < 0 || nm > max_nm || !valid_cigar($6)) {
        print $1 "\t*\t0\t*";
      } else {
        print $1 "\t" $3 "\t" $4 "\t" $6;
      }
    }
  ' "$input_path" > "$output_path"
}

compare_against_bwa() {
  local mapper_tsv="$1"
  local bwa_tsv="$2"
  local report_path="$3"

  awk '
    BEGIN { FS = "\t" }
    FNR == NR {
      mapper_chr[$1] = $2;
      mapper_pos[$1] = $3 + 0;
      mapper_cigar[$1] = $4;
      next;
    }
    {
      name = $1;
      bwa_chr = $2;
      bwa_pos = $3 + 0;
      bwa_cigar = $4;
      total++;

      mapper_is_mapped = ((name in mapper_chr) && mapper_chr[name] != "*");
      bwa_is_mapped = (bwa_chr != "*");

      if (mapper_is_mapped) mapper_mapped++;
      if (bwa_is_mapped) bwa_mapped++;

      if (mapper_is_mapped && bwa_is_mapped) {
        shared++;
        if (mapper_chr[name] == bwa_chr) {
          chrom_match++;
          delta = mapper_pos[name] - bwa_pos;
          if (delta < 0) delta = -delta;
          if (delta == 0) exact_pos++;
          if (delta <= 10) within10++;
          if (mapper_cigar[name] == bwa_cigar) exact_cigar++;
          if (delta == 0 && mapper_cigar[name] == bwa_cigar) exact_primary++;
        }
      } else if (mapper_is_mapped && !bwa_is_mapped) {
        mapper_only++;
      } else if (!mapper_is_mapped && bwa_is_mapped) {
        bwa_only++;
      }
    }
    END {
      printf "subset_reads=%d\n", total + 0;
      printf "mapper_mapped=%d\n", mapper_mapped + 0;
      printf "bwa_mapped=%d\n", bwa_mapped + 0;
      printf "shared_mapped=%d\n", shared + 0;
      printf "chrom_match=%d\n", chrom_match + 0;
      printf "within_10bp=%d\n", within10 + 0;
      printf "exact_pos=%d\n", exact_pos + 0;
      printf "exact_cigar=%d\n", exact_cigar + 0;
      printf "exact_primary=%d\n", exact_primary + 0;
      printf "mapper_only=%d\n", mapper_only + 0;
      printf "bwa_only=%d\n", bwa_only + 0;
    }
  ' "$mapper_tsv" "$bwa_tsv" > "$report_path"
}

make -C "$ROOT_DIR" -j4
RESOLVED_TIME_STYLE="$(detect_time_style)"

SUBSET="$OUT_DIR/reads_subset.fastq"
extract_subset "$READS" "$BENCH_READS" "$SUBSET"
READ_COUNT="$(count_reads "$SUBSET")"

MAPPER_STDOUT="$OUT_DIR/mapper.stdout.log"
MAPPER_STDERR="$OUT_DIR/mapper.stderr.log"
MAPPER_TIME="$OUT_DIR/mapper.time.log"
MAPPER_SAM="$OUT_DIR/mapper.sam"
MAPPER_STATS="$OUT_DIR/mapper.stats"
MAPPER_TSV="$OUT_DIR/mapper.tsv"

log "running mapper-speed on $READ_COUNT read(s)"
MAPPER_AVG_REAL="$(run_avg "$RUNS" "$ROOT_DIR/mapper" -I "$INDEX" -i "$SUBSET" -o "$OUT_DIR/mapper.avg.sam" -k "$K")"
run_once "mapper" "$MAPPER_STDOUT" "$MAPPER_STDERR" "$MAPPER_TIME" \
  "$ROOT_DIR/mapper" -I "$INDEX" -i "$SUBSET" -o "$MAPPER_SAM" -k "$K"
parse_time_report "$MAPPER_TIME" > "$OUT_DIR/mapper.metrics"
count_mapper_stats "$MAPPER_SAM" "$MAPPER_STATS"
convert_mapper_output "$MAPPER_SAM" "$MAPPER_TSV"

source "$OUT_DIR/mapper.metrics"
mapper_real_last="$real"
mapper_user="$user"
mapper_sys="$sys"
mapper_rss="$rss"
mapper_pf="$pf"
if [[ "$RESOLVED_TIME_STYLE" != "full" || "$mapper_rss" == "0" ]]; then
  mapper_rss="$(parse_mapper_peak_rss_bytes "$MAPPER_STDERR")"
fi
source "$MAPPER_STATS"
mapper_reads="$reads"
mapper_mapped="$mapped"
mapper_unmapped="$unmapped"
mapper_with_alt="$with_alt"
mapper_rps="$(awk -v reads="$READ_COUNT" -v real="$MAPPER_AVG_REAL" 'BEGIN { if (real > 0) printf "%.2f", reads / real; else print "0.00" }')"

BWA_AVAILABLE=0
if [[ "$COMPARE_BWA" == "1" ]]; then
  if command -v "$BWA_BIN" >/dev/null 2>&1; then
    BWA_AVAILABLE=1
  elif [[ -x "$BWA_BIN" ]]; then
    BWA_AVAILABLE=1
  fi
fi

if [[ "$BWA_AVAILABLE" == "1" ]]; then
  if [[ ! -f "${BWA_PREFIX}.bwt" ]]; then
    log "building bwa index"
    run_once "bwa-index" "$OUT_DIR/bwa.index.stdout.log" "$OUT_DIR/bwa.index.stderr.log" "$OUT_DIR/bwa.index.time.log" \
      "$BWA_BIN" index -p "$BWA_PREFIX" "$REF"
  fi

  BWA_STDOUT="$OUT_DIR/bwa.stdout.log"
  BWA_STDERR="$OUT_DIR/bwa.stderr.log"
  BWA_TIME="$OUT_DIR/bwa.time.log"
  BWA_SAM="$OUT_DIR/bwa.sam"
  BWA_STATS="$OUT_DIR/bwa.stats"
  BWA_TSV="$OUT_DIR/bwa.tsv"
  CORRECTNESS="$OUT_DIR/correctness.txt"

  log "running bwa mem on $READ_COUNT read(s)"
  BWA_AVG_REAL="$(run_avg "$RUNS" "$BWA_BIN" mem -t "$BWA_THREADS" "$BWA_PREFIX" "$SUBSET")"
  run_once "bwa" "$BWA_SAM" "$BWA_STDERR" "$BWA_TIME" \
    "$BWA_BIN" mem -t "$BWA_THREADS" "$BWA_PREFIX" "$SUBSET"
  parse_time_report "$BWA_TIME" > "$OUT_DIR/bwa.metrics"
  count_bwa_stats "$BWA_SAM" "$K" "$BWA_STATS"
  convert_bwa_output "$BWA_SAM" "$BWA_TSV" "$K"
  compare_against_bwa "$MAPPER_TSV" "$BWA_TSV" "$CORRECTNESS"

  source "$OUT_DIR/bwa.metrics"
  bwa_real_last="$real"
  bwa_user="$user"
  bwa_sys="$sys"
  bwa_rss="$rss"
  bwa_pf="$pf"
  source "$BWA_STATS"
  bwa_reads="$reads"
  bwa_mapped="$mapped"
  bwa_unmapped="$unmapped"
  bwa_with_alt="$with_alt"
  bwa_rps="$(awk -v reads="$READ_COUNT" -v real="$BWA_AVG_REAL" 'BEGIN { if (real > 0) printf "%.2f", reads / real; else print "0.00" }')"
  source "$CORRECTNESS"
fi

{
  echo "Benchmark Summary"
  echo "================="
  echo "reference: $REF"
  echo "reads: $READS"
  echo "bench_reads: $READ_COUNT"
  echo "index: $INDEX"
  echo "k: $K"
  echo "timing_backend: $RESOLVED_TIME_STYLE"
  echo
  echo "mapper-speed"
  echo "------------"
  echo "real_seconds_avg: $MAPPER_AVG_REAL"
  echo "real_seconds_last: $mapper_real_last"
  echo "user_seconds: $mapper_user"
  echo "sys_seconds: $mapper_sys"
  echo "reads_per_second: $mapper_rps"
  echo "peak_rss_bytes: $mapper_rss"
  echo "page_faults: $mapper_pf"
  echo "mapped: $mapper_mapped"
  echo "unmapped: $mapper_unmapped"
  echo "with_alt: $mapper_with_alt"

  if [[ "$BWA_AVAILABLE" == "1" ]]; then
    echo
    echo "bwa mem"
    echo "-------"
    echo "bwa_bin: $BWA_BIN"
    echo "real_seconds_avg: $BWA_AVG_REAL"
    echo "real_seconds_last: $bwa_real_last"
    echo "user_seconds: $bwa_user"
    echo "sys_seconds: $bwa_sys"
    echo "reads_per_second: $bwa_rps"
    echo "peak_rss_bytes: $bwa_rss"
    echo "page_faults: $bwa_pf"
    echo "mapped: $bwa_mapped"
    echo "unmapped: $bwa_unmapped"
    echo "with_alt: $bwa_with_alt"
    echo
    echo "correctness_vs_bwa"
    echo "------------------"
    echo "subset_reads: $subset_reads"
    echo "shared_mapped: $shared_mapped"
    echo "chrom_match: $chrom_match"
    echo "within_10bp: $within_10bp"
    echo "exact_pos: $exact_pos"
    echo "exact_cigar: $exact_cigar"
    echo "exact_primary: $exact_primary"
    echo "mapper_only: $mapper_only"
    echo "bwa_only: $bwa_only"
  else
    echo
    echo "bwa mem"
    echo "-------"
    echo "status: unavailable"
    echo "reason: bwa executable not found from this shell. Set BWA_BIN=/absolute/path/to/bwa and rerun."
  fi
} | tee "$OUT_DIR/summary.txt"

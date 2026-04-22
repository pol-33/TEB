#!/usr/bin/env bash
#
#SBATCH --job-name=mapper-speed-bench
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_DATA_DIR="${TMPDIR:-$ROOT_DIR}"

DENSE_INDEX="${DENSE_INDEX:-$DEFAULT_DATA_DIR/genome.idx}"
COMPACT_INDEX="${COMPACT_INDEX:-$DEFAULT_DATA_DIR/genome.compact.idx}"
READS="${READS:-$DEFAULT_DATA_DIR/reads_1M.fastq}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/bench-mn5-results}"
K="${K:-1}"
BUILD="${BUILD:-1}"
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
MAPPER_BIN="${MAPPER_BIN:-$ROOT_DIR/mapper}"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
BENCH_READS="${BENCH_READS:-0}"

mkdir -p "$OUT_DIR"

log() {
  printf '[bench-mn5] %s\n' "$*"
}

die() {
  printf '[bench-mn5][error] %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: ./bench-mn5.sh

This script runs four mapper-speed cases:
  1. dense index + native SIMD dispatch
  2. compact index + native SIMD dispatch
  3. dense index + AVX512 disabled
  4. compact index + AVX512 disabled

Environment overrides:
  DENSE_INDEX=...     Dense mapper index
  COMPACT_INDEX=...   Compact mapper index
  READS=...           Input FASTQ
  OUT_DIR=...         Output directory
  K=1                 Maximum edit distance
  BENCH_READS=0       Take the first N reads; 0 means all reads
  BUILD=1             Rebuild mapper before benchmarking
  BUILD_JOBS=4        make -j value
  MAPPER_BIN=...      Path to mapper executable
  TIME_BIN=/usr/bin/time
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ -f "$DENSE_INDEX" ]] || die "dense index not found: $DENSE_INDEX"
[[ -f "$COMPACT_INDEX" ]] || die "compact index not found: $COMPACT_INDEX"
[[ -f "$READS" ]] || die "reads FASTQ not found: $READS"
[[ "$K" =~ ^[0-9]+$ ]] || die "K must be an integer"
[[ "$BENCH_READS" =~ ^[0-9]+$ ]] || die "BENCH_READS must be an integer"

if [[ "$BUILD" == "1" ]]; then
  log "building mapper-speed"
  make -C "$ROOT_DIR" -j"$BUILD_JOBS"
fi

[[ -x "$MAPPER_BIN" ]] || die "mapper executable not found or not executable: $MAPPER_BIN"

detect_time_style() {
  local probe
  probe="$(mktemp)"
  if "$TIME_BIN" --version >/dev/null 2>&1; then
    rm -f "$probe"
    printf 'gnu\n'
    return
  fi
  if "$TIME_BIN" -l -o "$probe" true >/dev/null 2>&1; then
    rm -f "$probe"
    printf 'bsd\n'
    return
  fi
  rm -f "$probe"
  printf 'portable\n'
}

RESOLVED_TIME_STYLE="$(detect_time_style)"

count_reads() {
  awk 'END { printf "%d\n", NR / 4 }' "$1"
}

extract_subset() {
  local input_path="$1"
  local read_count="$2"
  local output_path="$3"
  if [[ "$read_count" -le 0 ]]; then
    ln -sf "$input_path" "$output_path"
  else
    awk -v limit="$((read_count * 4))" 'NR <= limit { print }' "$input_path" > "$output_path"
  fi
}

parse_time_report() {
  local time_path="$1"
  local style="$2"
  awk -v style="$style" '
    /^real=/                             { real = substr($0, 6) }
    /^user=/                             { user = substr($0, 6) }
    /^sys=/                              { sys = substr($0, 5) }
    /^rss_kb=/                           { rss = (substr($0, 8) + 0) * 1024 }
    /^pf=/                               { pf = substr($0, 4) + 0 }
    /real/ && /user/ && /sys/            { real = $1; user = $3; sys = $5 }
    /^real[[:space:]]+/                  { real = $2 }
    /^user[[:space:]]+/                  { user = $2 }
    /^sys[[:space:]]+/                   { sys = $2 }
    /maximum resident set size/          { rss = $1 }
    /page faults/                        { pf = $1 }
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

parse_active_simd() {
  local stderr_path="$1"
  awk '
    /single-thread verifier dispatch ready/ {
      line = $0;
      sub(/^.*\(simd=/, "", line);
      sub(/\).*$/, "", line);
      simd = line;
    }
    END {
      if (simd == "") simd = "unknown";
      print simd;
    }
  ' "$stderr_path"
}

run_once() {
  local stdout_path="$1"
  local stderr_path="$2"
  local time_path="$3"
  shift 3

  case "$RESOLVED_TIME_STYLE" in
    gnu)
      "$TIME_BIN" -f $'real=%e\nuser=%U\nsys=%S\nrss_kb=%M\npf=%F' -o "$time_path" \
        "$@" > "$stdout_path" 2> "$stderr_path"
      ;;
    bsd)
      "$TIME_BIN" -l -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
      ;;
    *)
      "$TIME_BIN" -p -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
      ;;
  esac
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

convert_mapper_output() {
  local input_path="$1"
  local output_path="$2"
  awk 'BEGIN { FS = "[ \t]+" } { print $1 "\t" $2 "\t" $3 "\t" $4 }' "$input_path" > "$output_path"
}

compare_mapper_outputs() {
  local lhs_tsv="$1"
  local rhs_tsv="$2"
  local report_path="$3"
  awk '
    BEGIN { FS = "\t" }
    FNR == NR {
      lhs_chr[$1] = $2;
      lhs_pos[$1] = $3 + 0;
      lhs_cigar[$1] = $4;
      next;
    }
    {
      name = $1;
      rhs_chr = $2;
      rhs_pos = $3 + 0;
      rhs_cigar = $4;
      total++;

      lhs_is_mapped = ((name in lhs_chr) && lhs_chr[name] != "*");
      rhs_is_mapped = (rhs_chr != "*");

      if (lhs_is_mapped) lhs_mapped++;
      if (rhs_is_mapped) rhs_mapped++;

      if (lhs_is_mapped && rhs_is_mapped) {
        shared++;
        if (lhs_chr[name] == rhs_chr) same_chrom++;
        if (lhs_chr[name] == rhs_chr && lhs_pos[name] == rhs_pos) same_pos++;
        if (lhs_chr[name] == rhs_chr && lhs_pos[name] == rhs_pos && lhs_cigar[name] == rhs_cigar) {
          exact_primary++;
        }
      } else if (lhs_is_mapped && !rhs_is_mapped) {
        lhs_only++;
      } else if (!lhs_is_mapped && rhs_is_mapped) {
        rhs_only++;
      }
    }
    END {
      printf "reads=%d\n", total + 0;
      printf "lhs_mapped=%d\n", lhs_mapped + 0;
      printf "rhs_mapped=%d\n", rhs_mapped + 0;
      printf "shared_mapped=%d\n", shared + 0;
      printf "same_chrom=%d\n", same_chrom + 0;
      printf "same_pos=%d\n", same_pos + 0;
      printf "exact_primary=%d\n", exact_primary + 0;
      printf "lhs_only_mapped=%d\n", lhs_only + 0;
      printf "rhs_only_mapped=%d\n", rhs_only + 0;
    }
  ' "$lhs_tsv" "$rhs_tsv" > "$report_path"
}

READS_SUBSET="$OUT_DIR/reads_subset.fastq"
extract_subset "$READS" "$BENCH_READS" "$READS_SUBSET"
READ_COUNT="$(count_reads "$READS_SUBSET")"

CASE_IDS=()
CASE_INDEX_KIND=()
CASE_SIMD_LABEL=()
CASE_SIMD_REQUEST=()
CASE_SIMD_EFFECTIVE=()
CASE_REAL=()
CASE_USER=()
CASE_SYS=()
CASE_RSS=()
CASE_MAPPED=()
CASE_UNMAPPED=()
CASE_WITH_ALT=()
CASE_DIRS=()

run_case() {
  local case_id="$1"
  local index_kind="$2"
  local index_path="$3"
  local simd_label="$4"
  local simd_request="$5"
  local case_dir="$OUT_DIR/$case_id"
  local stdout_path="$case_dir/mapper.stdout.log"
  local stderr_path="$case_dir/mapper.stderr.log"
  local time_path="$case_dir/mapper.time.log"
  local sam_path="$case_dir/mapper.sam"
  local stats_path="$case_dir/mapper.stats"
  local tsv_path="$case_dir/mapper.tsv"
  local metrics_path="$case_dir/mapper.metrics"

  mkdir -p "$case_dir"
  log "running $case_id"

  if [[ "$simd_request" == "native" ]]; then
    run_once "$stdout_path" "$stderr_path" "$time_path" \
      "$MAPPER_BIN" -I "$index_path" -i "$READS_SUBSET" -o "$sam_path" -k "$K"
  else
    run_once "$stdout_path" "$stderr_path" "$time_path" \
      env MAPPER_SPEED_MAX_SIMD="$simd_request" "$MAPPER_BIN" \
      -I "$index_path" -i "$READS_SUBSET" -o "$sam_path" -k "$K"
  fi

  parse_time_report "$time_path" "$RESOLVED_TIME_STYLE" > "$metrics_path"
  count_mapper_stats "$sam_path" "$stats_path"
  convert_mapper_output "$sam_path" "$tsv_path"

  source "$metrics_path"
  local real_seconds="$real"
  local user_seconds="$user"
  local sys_seconds="$sys"
  local peak_rss_bytes="$rss"

  if [[ "$peak_rss_bytes" == "0" ]]; then
    peak_rss_bytes="$(parse_mapper_peak_rss_bytes "$stderr_path")"
  fi

  source "$stats_path"
  local effective_simd
  effective_simd="$(parse_active_simd "$stderr_path")"

  CASE_IDS+=("$case_id")
  CASE_INDEX_KIND+=("$index_kind")
  CASE_SIMD_LABEL+=("$simd_label")
  CASE_SIMD_REQUEST+=("$simd_request")
  CASE_SIMD_EFFECTIVE+=("$effective_simd")
  CASE_REAL+=("$real_seconds")
  CASE_USER+=("$user_seconds")
  CASE_SYS+=("$sys_seconds")
  CASE_RSS+=("$peak_rss_bytes")
  CASE_MAPPED+=("$mapped")
  CASE_UNMAPPED+=("$unmapped")
  CASE_WITH_ALT+=("$with_alt")
  CASE_DIRS+=("$case_dir")

  {
    echo "case_id=$case_id"
    echo "index_kind=$index_kind"
    echo "simd_label=$simd_label"
    echo "simd_request=$simd_request"
    echo "effective_simd=$effective_simd"
    echo "real_seconds=$real_seconds"
    echo "user_seconds=$user_seconds"
    echo "sys_seconds=$sys_seconds"
    echo "peak_rss_bytes=$peak_rss_bytes"
    echo "mapped=$mapped"
    echo "unmapped=$unmapped"
    echo "with_alt=$with_alt"
  } > "$case_dir/run.meta"
}

run_case "dense-native" "dense" "$DENSE_INDEX" "native" "native"
run_case "compact-native" "compact" "$COMPACT_INDEX" "native" "native"
run_case "dense-noavx512" "dense" "$DENSE_INDEX" "avx512_off" "avx2"
run_case "compact-noavx512" "compact" "$COMPACT_INDEX" "avx512_off" "avx2"

RUNS_TSV="$OUT_DIR/runs.tsv"
{
  printf "case_id\tindex_kind\tsimd_label\tsimd_request\teffective_simd\treal_seconds\tuser_seconds\tsys_seconds\tpeak_rss_bytes\tmapped\tunmapped\twith_alt\tcase_dir\n"
  for i in "${!CASE_IDS[@]}"; do
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${CASE_IDS[$i]}" "${CASE_INDEX_KIND[$i]}" "${CASE_SIMD_LABEL[$i]}" \
      "${CASE_SIMD_REQUEST[$i]}" "${CASE_SIMD_EFFECTIVE[$i]}" \
      "${CASE_REAL[$i]}" "${CASE_USER[$i]}" "${CASE_SYS[$i]}" \
      "${CASE_RSS[$i]}" "${CASE_MAPPED[$i]}" "${CASE_UNMAPPED[$i]}" \
      "${CASE_WITH_ALT[$i]}" "${CASE_DIRS[$i]}"
  done
} > "$RUNS_TSV"

COMPARISONS_TSV="$OUT_DIR/comparisons.tsv"
printf "lhs_case\trhs_case\tshared_mapped\texact_primary\tlhs_only_mapped\trhs_only_mapped\tdelta_mapped\tdelta_unmapped\tdelta_with_alt\tsame_chrom\tsame_pos\tcomparison_file\n" > "$COMPARISONS_TSV"

for ((i = 0; i < ${#CASE_IDS[@]}; ++i)); do
  for ((j = i + 1; j < ${#CASE_IDS[@]}; ++j)); do
    lhs_case="${CASE_IDS[$i]}"
    rhs_case="${CASE_IDS[$j]}"
    comparison_file="$OUT_DIR/compare-${lhs_case}-vs-${rhs_case}.txt"
    compare_mapper_outputs "${CASE_DIRS[$i]}/mapper.tsv" "${CASE_DIRS[$j]}/mapper.tsv" "$comparison_file"
    source "$comparison_file"
    delta_mapped="$(awk -v a="${CASE_MAPPED[$i]}" -v b="${CASE_MAPPED[$j]}" 'BEGIN { print a - b }')"
    delta_unmapped="$(awk -v a="${CASE_UNMAPPED[$i]}" -v b="${CASE_UNMAPPED[$j]}" 'BEGIN { print a - b }')"
    delta_with_alt="$(awk -v a="${CASE_WITH_ALT[$i]}" -v b="${CASE_WITH_ALT[$j]}" 'BEGIN { print a - b }')"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "$lhs_case" "$rhs_case" "$shared_mapped" "$exact_primary" \
      "$lhs_only_mapped" "$rhs_only_mapped" "$delta_mapped" \
      "$delta_unmapped" "$delta_with_alt" "$same_chrom" "$same_pos" \
      "$comparison_file" >> "$COMPARISONS_TSV"
  done
done

SUMMARY_TXT="$OUT_DIR/summary.txt"
{
  echo "MN5 Mapper-Speed Benchmark"
  echo "========================="
  echo "date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "host: $(hostname)"
  echo "slurm_job_id: ${SLURM_JOB_ID:-n/a}"
  echo "reads: $READS"
  echo "bench_reads: $READ_COUNT"
  echo "dense_index: $DENSE_INDEX"
  echo "compact_index: $COMPACT_INDEX"
  echo "k: $K"
  echo "timing_backend: $RESOLVED_TIME_STYLE"
  echo
  echo "Runs"
  echo "----"
  printf "%-20s %-8s %-12s %-20s %12s %16s %10s %10s %10s\n" \
    "case" "index" "simd" "effective" "real_s" "peak_rss_bytes" "mapped" "unmapped" "with_alt"
  for i in "${!CASE_IDS[@]}"; do
    printf "%-20s %-8s %-12s %-20s %12s %16s %10s %10s %10s\n" \
      "${CASE_IDS[$i]}" "${CASE_INDEX_KIND[$i]}" "${CASE_SIMD_LABEL[$i]}" \
      "${CASE_SIMD_EFFECTIVE[$i]}" "${CASE_REAL[$i]}" "${CASE_RSS[$i]}" \
      "${CASE_MAPPED[$i]}" "${CASE_UNMAPPED[$i]}" "${CASE_WITH_ALT[$i]}"
  done
  echo
  echo "Pairwise Comparisons"
  echo "--------------------"
  printf "%-20s %-20s %14s %14s %14s %14s %14s %14s %14s\n" \
    "lhs" "rhs" "shared" "exact" "lhs_only" "rhs_only" "d_mapped" "d_unmapped" "d_with_alt"
  while IFS=$'\t' read -r lhs_case rhs_case shared_mapped exact_primary lhs_only_mapped rhs_only_mapped delta_mapped delta_unmapped delta_with_alt same_chrom same_pos comparison_file; do
    [[ "$lhs_case" == "lhs_case" ]] && continue
    printf "%-20s %-20s %14s %14s %14s %14s %14s %14s %14s\n" \
      "$lhs_case" "$rhs_case" "$shared_mapped" "$exact_primary" "$lhs_only_mapped" \
      "$rhs_only_mapped" "$delta_mapped" "$delta_unmapped" "$delta_with_alt"
  done < "$COMPARISONS_TSV"
  echo
  echo "machine_readable_runs: $RUNS_TSV"
  echo "machine_readable_comparisons: $COMPARISONS_TSV"
} | tee "$SUMMARY_TXT"

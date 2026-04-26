#!/usr/bin/env bash
#
#SBATCH --job-name=mapper-speed-bench
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --account=nct_370
#SBATCH --qos=acc_debug

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_DATA_DIR="${TMPDIR:-$ROOT_DIR}"

REF="${REF:-$DEFAULT_DATA_DIR/genome.fa}"
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
COMPARE_BWA="${COMPARE_BWA:-1}"
BWA_MODULE="${BWA_MODULE:-bwa/0.7.17}"
BWA_BIN="${BWA_BIN:-bwa}"
BWA_THREADS="${BWA_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
BWA_PREFIX="${BWA_PREFIX:-$OUT_DIR/reference.bwa}"

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
  REF=...             Reference FASTA for bwa index
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
  COMPARE_BWA=1       Also run bwa mem and compare against it
  BWA_MODULE=...      Module name to load when bwa is not already in PATH
  BWA_BIN=bwa         bwa executable name or absolute path
  BWA_THREADS=1       Thread count for bwa mem
  BWA_PREFIX=...      Prefix for bwa index files
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
[[ "$COMPARE_BWA" =~ ^[01]$ ]] || die "COMPARE_BWA must be 0 or 1"
[[ "$BWA_THREADS" =~ ^[0-9]+$ ]] || die "BWA_THREADS must be an integer"
if [[ "$COMPARE_BWA" == "1" ]]; then
  [[ -f "$REF" ]] || die "reference FASTA not found: $REF"
fi

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
    local input_abs
    input_abs="$(cd "$(dirname "$input_path")" && pwd)/$(basename "$input_path")"
    ln -sf "$input_abs" "$output_path"
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

absolute_path() {
  local path="$1"
  if [[ "$path" = /* ]]; then
    printf '%s\n' "$path"
  else
    printf '%s/%s\n' "$(cd "$(dirname "$path")" && pwd)" "$(basename "$path")"
  fi
}

ensure_module_cmd() {
  if command -v module >/dev/null 2>&1; then
    return 0
  fi
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! command -v module >/dev/null 2>&1 && [[ -f /usr/share/lmod/lmod/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/lmod/lmod/init/bash
  fi
  if ! command -v module >/dev/null 2>&1 && [[ -f /apps/modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /apps/modules/init/bash
  fi
  command -v module >/dev/null 2>&1
}

bwa_index_complete() {
  local prefix="$1"
  local ext
  for ext in amb ann bwt pac sa; do
    [[ -f "${prefix}.${ext}" ]] || return 1
  done
  return 0
}

log_bwa_index_state() {
  local prefix="$1"
  local ext
  local found=()
  local missing=()
  for ext in amb ann bwt pac sa; do
    if [[ -f "${prefix}.${ext}" ]]; then
      found+=("${prefix}.${ext}")
    else
      missing+=("${prefix}.${ext}")
    fi
  done
  if [[ "${#found[@]}" -gt 0 ]]; then
    log "existing bwa index parts: ${found[*]}"
  fi
  if [[ "${#missing[@]}" -gt 0 ]]; then
    log "missing bwa index parts: ${missing[*]}"
  fi
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
      printf "reads=%d\n", total + 0;
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

mapper_case_complete() {
  local case_dir="$1"
  [[ -f "$case_dir/mapper.sam" &&
     -f "$case_dir/mapper.stats" &&
     -f "$case_dir/mapper.tsv" &&
     -f "$case_dir/mapper.metrics" &&
     -f "$case_dir/mapper.stderr.log" ]]
}

bwa_run_complete() {
  [[ -f "$OUT_DIR/bwa.sam" &&
     -f "$OUT_DIR/bwa.stats" &&
     -f "$OUT_DIR/bwa.tsv" &&
     -f "$OUT_DIR/bwa.metrics" &&
     -f "$OUT_DIR/bwa.stderr.log" ]]
}

load_case_results() {
  local case_id="$1"
  local index_kind="$2"
  local simd_label="$3"
  local simd_request="$4"
  local case_dir="$OUT_DIR/$case_id"
  local stderr_path="$case_dir/mapper.stderr.log"

  source "$case_dir/mapper.metrics"
  local real_seconds="$real"
  local user_seconds="$user"
  local sys_seconds="$sys"
  local peak_rss_bytes="$rss"

  if [[ "$peak_rss_bytes" == "0" ]]; then
    peak_rss_bytes="$(parse_mapper_peak_rss_bytes "$stderr_path")"
  fi

  source "$case_dir/mapper.stats"
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
  if mapper_case_complete "$case_dir"; then
    log "skipping $case_id (reusing existing results)"
    load_case_results "$case_id" "$index_kind" "$simd_label" "$simd_request"
    return 0
  fi
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
  load_case_results "$case_id" "$index_kind" "$simd_label" "$simd_request"
}

run_case "dense-native" "dense" "$DENSE_INDEX" "native" "native"
run_case "compact-native" "compact" "$COMPACT_INDEX" "native" "native"
run_case "dense-noavx512" "dense" "$DENSE_INDEX" "avx512_off" "avx2"
run_case "compact-noavx512" "compact" "$COMPACT_INDEX" "avx512_off" "avx2"

BWA_AVAILABLE=0
BWA_LOADED_MODULE="n/a"
bwa_real_seconds="n/a"
bwa_user_seconds="n/a"
bwa_sys_seconds="n/a"
bwa_peak_rss_bytes="n/a"
bwa_mapped="n/a"
bwa_unmapped="n/a"
bwa_with_alt="n/a"

ensure_bwa_available() {
  if command -v "$BWA_BIN" >/dev/null 2>&1; then
    BWA_BIN="$(command -v "$BWA_BIN")"
    BWA_AVAILABLE=1
    BWA_LOADED_MODULE="already_in_path"
    return 0
  fi
  if [[ -x "$BWA_BIN" ]]; then
    BWA_BIN="$(absolute_path "$BWA_BIN")"
    BWA_AVAILABLE=1
    BWA_LOADED_MODULE="custom_path"
    return 0
  fi
  if ! ensure_module_cmd; then
    return 1
  fi
  log "loading module $BWA_MODULE"
  module load "$BWA_MODULE"
  hash -r
  if command -v "$BWA_BIN" >/dev/null 2>&1; then
    BWA_BIN="$(command -v "$BWA_BIN")"
    BWA_AVAILABLE=1
    BWA_LOADED_MODULE="$BWA_MODULE"
    return 0
  fi
  return 1
}

if [[ "$COMPARE_BWA" == "1" ]]; then
  ensure_bwa_available || die "bwa is unavailable. Load it manually or set BWA_MODULE/BWA_BIN correctly."

  if ! bwa_index_complete "$BWA_PREFIX"; then
    log_bwa_index_state "$BWA_PREFIX"
    log "building bwa index (the final .sa stage can take a long time on full GRCh38)"
    run_once "$OUT_DIR/bwa.index.stdout.log" "$OUT_DIR/bwa.index.stderr.log" "$OUT_DIR/bwa.index.time.log" \
      "$BWA_BIN" index -p "$BWA_PREFIX" "$REF"
    bwa_index_complete "$BWA_PREFIX" || die "bwa index did not finish correctly for prefix: $BWA_PREFIX"
  else
    log "reusing existing bwa index at $BWA_PREFIX"
  fi

  BWA_STDERR="$OUT_DIR/bwa.stderr.log"
  BWA_TIME="$OUT_DIR/bwa.time.log"
  BWA_SAM="$OUT_DIR/bwa.sam"
  BWA_STATS="$OUT_DIR/bwa.stats"
  BWA_TSV="$OUT_DIR/bwa.tsv"
  BWA_METRICS="$OUT_DIR/bwa.metrics"

  if bwa_run_complete; then
    log "skipping bwa mem (reusing existing results)"
  else
    log "running bwa mem with $BWA_BIN"
    run_once "$BWA_SAM" "$BWA_STDERR" "$BWA_TIME" \
      "$BWA_BIN" mem -t "$BWA_THREADS" "$BWA_PREFIX" "$READS_SUBSET"
    parse_time_report "$BWA_TIME" "$RESOLVED_TIME_STYLE" > "$BWA_METRICS"
    count_bwa_stats "$BWA_SAM" "$K" "$BWA_STATS"
    convert_bwa_output "$BWA_SAM" "$BWA_TSV" "$K"
  fi

  source "$BWA_METRICS"
  bwa_real_seconds="$real"
  bwa_user_seconds="$user"
  bwa_sys_seconds="$sys"
  bwa_peak_rss_bytes="$rss"

  source "$BWA_STATS"
  bwa_mapped="$mapped"
  bwa_unmapped="$unmapped"
  bwa_with_alt="$with_alt"
fi

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
  if [[ "$BWA_AVAILABLE" == "1" ]]; then
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "bwa-mem" "baseline" "baseline" "bwa-mem" "bwa-mem" \
      "$bwa_real_seconds" "$bwa_user_seconds" "$bwa_sys_seconds" \
      "$bwa_peak_rss_bytes" "$bwa_mapped" "$bwa_unmapped" \
      "$bwa_with_alt" "$OUT_DIR"
  fi
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

BWA_CORRECTNESS_TSV="$OUT_DIR/bwa-correctness.tsv"
if [[ "$BWA_AVAILABLE" == "1" ]]; then
  printf "mapper_case\tshared_mapped\tchrom_match\twithin_10bp\texact_pos\texact_cigar\texact_primary\tmapper_only\tbwa_only\tdelta_mapped\tdelta_unmapped\tdelta_with_alt\tcomparison_file\n" > "$BWA_CORRECTNESS_TSV"
  for i in "${!CASE_IDS[@]}"; do
    comparison_file="$OUT_DIR/correctness-${CASE_IDS[$i]}-vs-bwa.txt"
    compare_against_bwa "${CASE_DIRS[$i]}/mapper.tsv" "$BWA_TSV" "$comparison_file"
    source "$comparison_file"
    delta_mapped="$(awk -v a="${CASE_MAPPED[$i]}" -v b="$bwa_mapped" 'BEGIN { print a - b }')"
    delta_unmapped="$(awk -v a="${CASE_UNMAPPED[$i]}" -v b="$bwa_unmapped" 'BEGIN { print a - b }')"
    delta_with_alt="$(awk -v a="${CASE_WITH_ALT[$i]}" -v b="$bwa_with_alt" 'BEGIN { print a - b }')"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${CASE_IDS[$i]}" "$shared_mapped" "$chrom_match" "$within_10bp" \
      "$exact_pos" "$exact_cigar" "$exact_primary" "$mapper_only" \
      "$bwa_only" "$delta_mapped" "$delta_unmapped" "$delta_with_alt" \
      "$comparison_file" >> "$BWA_CORRECTNESS_TSV"
  done
fi

SUMMARY_TXT="$OUT_DIR/summary.txt"
{
  echo "MN5 Mapper-Speed Benchmark"
  echo "========================="
  echo "date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
  echo "host: $(hostname)"
  echo "slurm_job_id: ${SLURM_JOB_ID:-n/a}"
  echo "reference: $REF"
  echo "reads: $READS"
  echo "bench_reads: $READ_COUNT"
  echo "dense_index: $DENSE_INDEX"
  echo "compact_index: $COMPACT_INDEX"
  echo "k: $K"
  echo "timing_backend: $RESOLVED_TIME_STYLE"
  echo "bwa_module: $BWA_MODULE"
  echo "bwa_loaded_module: $BWA_LOADED_MODULE"
  echo "bwa_prefix: $BWA_PREFIX"
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
  if [[ "$BWA_AVAILABLE" == "1" ]]; then
    printf "%-20s %-8s %-12s %-20s %12s %16s %10s %10s %10s\n" \
      "bwa-mem" "bwa" "baseline" "bwa-mem" "$bwa_real_seconds" \
      "$bwa_peak_rss_bytes" "$bwa_mapped" "$bwa_unmapped" "$bwa_with_alt"
  fi
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
  if [[ "$BWA_AVAILABLE" == "1" ]]; then
    echo
    echo "Correctness Vs BWA"
    echo "------------------"
    printf "%-20s %12s %12s %12s %12s %12s %12s %12s %12s\n" \
      "mapper_case" "shared" "chrom" "within10" "exact_pos" "exact_cigar" "exact" "m_only" "bwa_only"
    while IFS=$'\t' read -r mapper_case shared_mapped chrom_match within_10bp exact_pos exact_cigar exact_primary mapper_only bwa_only delta_mapped delta_unmapped delta_with_alt comparison_file; do
      [[ "$mapper_case" == "mapper_case" ]] && continue
      printf "%-20s %12s %12s %12s %12s %12s %12s %12s %12s\n" \
        "$mapper_case" "$shared_mapped" "$chrom_match" "$within_10bp" \
        "$exact_pos" "$exact_cigar" "$exact_primary" "$mapper_only" "$bwa_only"
    done < "$BWA_CORRECTNESS_TSV"
  fi
  echo
  echo "machine_readable_runs: $RUNS_TSV"
  echo "machine_readable_comparisons: $COMPARISONS_TSV"
  if [[ "$BWA_AVAILABLE" == "1" ]]; then
    echo "machine_readable_bwa_correctness: $BWA_CORRECTNESS_TSV"
  fi
} | tee "$SUMMARY_TXT"

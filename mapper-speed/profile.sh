#!/usr/bin/env bash
#
#SBATCH --job-name=mapper-speed-profile
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

INDEX="${INDEX:-$DEFAULT_DATA_DIR/genome.compact.idx}"
READS="${READS:-$DEFAULT_DATA_DIR/reads_1M.fastq}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/profile-results}"
K="${K:-1}"
BUILD="${BUILD:-1}"
BUILD_JOBS="${BUILD_JOBS:-${SLURM_CPUS_PER_TASK:-4}}"
MAPPER_BIN="${MAPPER_BIN:-$ROOT_DIR/mapper}"
PROFILE_MODE="${PROFILE_MODE:-all}"
PROFILE_READS="${PROFILE_READS:-5000}"
PERF_BIN="${PERF_BIN:-perf}"
VALGRIND_BIN="${VALGRIND_BIN:-valgrind}"
CALLGRIND_ANNOTATE_BIN="${CALLGRIND_ANNOTATE_BIN:-callgrind_annotate}"
MS_PRINT_BIN="${MS_PRINT_BIN:-ms_print}"
PERF_STAT_EVENTS="${PERF_STAT_EVENTS:-task-clock,cycles,instructions,branches,branch-misses,cache-references,cache-misses,minor-faults,major-faults}"
PERF_RECORD_FREQ="${PERF_RECORD_FREQ:-999}"
PERF_REPORT_LINES="${PERF_REPORT_LINES:-25}"
CALLGRIND_ARGS="${CALLGRIND_ARGS:---dump-instr=no --collect-jumps=no --collect-systime=no}"
MASSIF_ARGS="${MASSIF_ARGS:---time-unit=i --stacks=no}"
VALGRIND_SAFE_BUILD="${VALGRIND_SAFE_BUILD:-1}"
VALGRIND_SIMD_CAP="${VALGRIND_SIMD_CAP:-generic}"
VALGRIND_MAPPER_BIN=""

mkdir -p "$OUT_DIR"

log() {
  printf '[profile] %s\n' "$*"
}

die() {
  printf '[profile][error] %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: ./profile.sh

This script profiles mapper-speed and writes:
  - raw perf / valgrind artifacts
  - a compact markdown summary that is easy to read or paste into an LLM

Environment overrides:
  INDEX=...                 Mapper index to profile
  READS=...                 Input FASTQ
  OUT_DIR=...               Output directory
  K=1                       Maximum edit distance
  PROFILE_READS=5000        Use the first N reads; 0 means all reads
  PROFILE_MODE=all          all|perf|valgrind|perf-stat|perf-record|callgrind|massif
  BUILD=1                   Rebuild mapper before profiling
  BUILD_JOBS=4              make -j value
  MAPPER_BIN=...            Path to mapper executable
  PERF_BIN=perf             perf executable
  VALGRIND_BIN=valgrind     valgrind executable
  PERF_STAT_EVENTS=...      perf stat events
  PERF_RECORD_FREQ=999      perf record sample frequency
  PERF_REPORT_LINES=25      Hotspot lines kept in summary
  CALLGRIND_ARGS=...        Extra callgrind arguments
  MASSIF_ARGS=...           Extra massif arguments
  VALGRIND_SAFE_BUILD=1     Rebuild a portable mapper for valgrind tools
  VALGRIND_SIMD_CAP=generic SIMD cap used under valgrind

Examples:
  PROFILE_MODE=perf PROFILE_READS=20000 ./profile.sh
  PROFILE_MODE=callgrind PROFILE_READS=1000 MAPPER_SPEED_MAX_SIMD=avx2 ./profile.sh
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ -f "$INDEX" ]] || die "index not found: $INDEX"
[[ -f "$READS" ]] || die "reads FASTQ not found: $READS"
[[ "$K" =~ ^[0-9]+$ ]] || die "K must be an integer"
[[ "$PROFILE_READS" =~ ^[0-9]+$ ]] || die "PROFILE_READS must be an integer"
[[ -x "$MAPPER_BIN" ]] || die "mapper executable not found or not executable: $MAPPER_BIN"

if [[ "$BUILD" == "1" ]]; then
  log "building mapper-speed"
  make -C "$ROOT_DIR" -j"$BUILD_JOBS"
fi

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

count_reads() {
  awk 'END { printf "%d\n", NR / 4 }' "$1"
}

mode_enabled() {
  local probe="$1"
  case "$PROFILE_MODE" in
    all) return 0 ;;
    perf) [[ "$probe" == perf-stat || "$probe" == perf-record ]] ;;
    valgrind) [[ "$probe" == callgrind || "$probe" == massif ]] ;;
    *) [[ "$PROFILE_MODE" == "$probe" ]] ;;
  esac
}

format_bytes_iec() {
  awk -v bytes="$1" '
    function abs(v) { return v < 0 ? -v : v }
    BEGIN {
      split("B KiB MiB GiB TiB", units, " ");
      value = bytes + 0.0;
      unit = 1;
      while (abs(value) >= 1024.0 && unit < 5) {
        value /= 1024.0;
        ++unit;
      }
      printf "%.2f %s", value, units[unit];
    }'
}

write_metadata() {
  local path="$1"
  {
    echo "date=$(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "host=$(hostname)"
    echo "slurm_job_id=${SLURM_JOB_ID:-n/a}"
    echo "index=$INDEX"
    echo "reads=$READS"
    echo "profile_reads=$PROFILE_READ_COUNT"
    echo "k=$K"
    echo "profile_mode=$PROFILE_MODE"
    echo "mapper_bin=$MAPPER_BIN"
    echo "simd_cap=${MAPPER_SPEED_MAX_SIMD:-native}"
    echo "disable_avx512=${MAPPER_SPEED_DISABLE_AVX512:-0}"
    echo "perf_bin=$PERF_BIN"
    echo "valgrind_bin=$VALGRIND_BIN"
  } > "$path"
}

append_markdown_header() {
  local summary_path="$1"
  {
    echo "# Mapper-Speed Profile Summary"
    echo
    echo "## Run"
    echo
    echo "| Field | Value |"
    echo "| --- | --- |"
    echo "| Date | $(date '+%Y-%m-%d %H:%M:%S %Z') |"
    echo "| Host | $(hostname) |"
    echo "| Slurm job id | ${SLURM_JOB_ID:-n/a} |"
    echo "| Index | \`$INDEX\` |"
    echo "| Reads | \`$READS\` |"
    echo "| Profile reads | $PROFILE_READ_COUNT |"
    echo "| k | $K |"
    echo "| Mode | \`$PROFILE_MODE\` |"
    echo "| SIMD cap | \`${MAPPER_SPEED_MAX_SIMD:-native}\` |"
    echo "| Disable AVX512 | \`${MAPPER_SPEED_DISABLE_AVX512:-0}\` |"
    echo
  } > "$summary_path"
}

append_tool_unavailable() {
  local summary_path="$1"
  local tool_name="$2"
  local detail="$3"
  {
    echo "## $tool_name"
    echo
    echo "$detail"
    echo
  } >> "$summary_path"
}

append_fenced_excerpt() {
  local summary_path="$1"
  local title="$2"
  local file_path="$3"
  local line_count="$4"
  {
    echo "### $title"
    echo
    echo '```text'
    sed -n "1,${line_count}p" "$file_path"
    echo '```'
    echo
  } >> "$summary_path"
}

perf_stat_section() {
  local raw_path="$1"
  local summary_path="$2"
  local markdown_path="$OUT_DIR/perf-stat/perf-stat-table.md"
  awk -F';' '
    function clean(v) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", v);
      return v;
    }
    function numeric(v) {
      gsub(/,/, "", v);
      return v;
    }
    BEGIN {
      print "| Event | Value | Unit |";
      print "| --- | --- | --- |";
    }
    NF >= 3 {
      value = clean($1);
      unit = clean($2);
      event = clean($3);
      if (event == "" || value == "" || value ~ /^</) {
        next;
      }
      metric[event] = numeric(value);
      metric_text[event] = value;
      unit_text[event] = (unit == "" ? "-" : unit);
      printf "| %s | %s | %s |\n", event, value, unit_text[event];
    }
    END {
      ipc = "-";
      branch_rate = "-";
      cache_rate = "-";
      if (metric["cycles"] > 0 && metric["instructions"] > 0) {
        ipc = sprintf("%.3f", metric["instructions"] / metric["cycles"]);
      }
      if (metric["branches"] > 0 && metric["branch-misses"] >= 0) {
        branch_rate = sprintf("%.2f%%", (100.0 * metric["branch-misses"]) / metric["branches"]);
      }
      if (metric["cache-references"] > 0 && metric["cache-misses"] >= 0) {
        cache_rate = sprintf("%.2f%%", (100.0 * metric["cache-misses"]) / metric["cache-references"]);
      }
      print "" > "/dev/stderr";
      print "ipc=" ipc > "/dev/stderr";
      print "branch_miss_rate=" branch_rate > "/dev/stderr";
      print "cache_miss_rate=" cache_rate > "/dev/stderr";
    }' "$raw_path" > "$markdown_path" 2> "$OUT_DIR/perf-stat/perf-stat-derived.env"

  # shellcheck disable=SC1090
  source "$OUT_DIR/perf-stat/perf-stat-derived.env"
  {
    echo "## Perf Stat"
    echo
    echo "| Derived metric | Value |"
    echo "| --- | --- |"
    echo "| IPC | ${ipc:-n/a} |"
    echo "| Branch miss rate | ${branch_miss_rate:-n/a} |"
    echo "| Cache miss rate | ${cache_miss_rate:-n/a} |"
    echo
    cat "$markdown_path"
    echo
  } >> "$summary_path"
}

perf_record_section() {
  local report_path="$1"
  local summary_path="$2"
  local top_path="$OUT_DIR/perf-record/perf-top.txt"
  awk -v limit="$PERF_REPORT_LINES" '
    BEGIN { count = 0 }
    /^[[:space:]]*[0-9]+\.[0-9]+%/ {
      print;
      ++count;
      if (count >= limit) exit;
    }' "$report_path" > "$top_path"
  {
    echo "## Perf Record"
    echo
    echo "Raw files:"
    echo
    echo "- \`$OUT_DIR/perf-record/perf.data\`"
    echo "- \`$OUT_DIR/perf-record/perf.report.txt\`"
    echo
  } >> "$summary_path"
  if [[ -s "$top_path" ]]; then
    append_fenced_excerpt "$summary_path" "Top Hotspots" "$top_path" "$PERF_REPORT_LINES"
  else
    echo "No hotspot lines were extracted from perf report." >> "$summary_path"
    echo >> "$summary_path"
  fi
}

callgrind_section() {
  local annotate_path="$1"
  local summary_path="$2"
  local top_path="$OUT_DIR/callgrind/callgrind-top.txt"
  awk -v limit="$PERF_REPORT_LINES" '
    /^[[:space:]]*[0-9,]+[[:space:]]+/ {
      print;
      ++count;
      if (count >= limit) exit;
    }' "$annotate_path" > "$top_path"
  {
    echo "## Callgrind"
    echo
    echo "Raw files:"
    echo
    echo "- \`$OUT_DIR/callgrind/callgrind.out\`"
    echo "- \`$OUT_DIR/callgrind/callgrind.annotate.txt\`"
    echo
  } >> "$summary_path"
  if [[ -s "$top_path" ]]; then
    append_fenced_excerpt "$summary_path" "Top Functions" "$top_path" "$PERF_REPORT_LINES"
  else
    echo "No function lines were extracted from callgrind output." >> "$summary_path"
    echo >> "$summary_path"
  fi
}

massif_peak_env() {
  local massif_path="$1"
  awk -F'=' '
    /^snapshot=/ { cur = $2 + 0; heap = 0; extra = 0; stacks = 0; next }
    /^mem_heap_B=/ { heap = $2 + 0; next }
    /^mem_heap_extra_B=/ { extra = $2 + 0; next }
    /^mem_stacks_B=/ {
      stacks = $2 + 0;
      total = heap + extra + stacks;
      if (total > best) {
        best = total;
        best_snapshot = cur;
        best_heap = heap;
        best_extra = extra;
        best_stacks = stacks;
      }
    }
    END {
      printf "peak_snapshot=%d\n", best_snapshot + 0;
      printf "peak_total_bytes=%d\n", best + 0;
      printf "peak_heap_bytes=%d\n", best_heap + 0;
      printf "peak_heap_extra_bytes=%d\n", best_extra + 0;
      printf "peak_stacks_bytes=%d\n", best_stacks + 0;
    }' "$massif_path"
}

massif_section() {
  local massif_path="$1"
  local summary_path="$2"
  local env_path="$OUT_DIR/massif/massif-peak.env"
  massif_peak_env "$massif_path" > "$env_path"
  # shellcheck disable=SC1090
  source "$env_path"
  {
    echo "## Massif"
    echo
    echo "| Metric | Value |"
    echo "| --- | --- |"
    echo "| Peak snapshot | ${peak_snapshot:-0} |"
    echo "| Peak total | $(format_bytes_iec "${peak_total_bytes:-0}") |"
    echo "| Peak heap | $(format_bytes_iec "${peak_heap_bytes:-0}") |"
    echo "| Peak heap extra | $(format_bytes_iec "${peak_heap_extra_bytes:-0}") |"
    echo "| Peak stacks | $(format_bytes_iec "${peak_stacks_bytes:-0}") |"
    echo
    echo "Raw files:"
    echo
    echo "- \`$OUT_DIR/massif/massif.out\`"
    if [[ -f "$OUT_DIR/massif/massif.ms_print.txt" ]]; then
      echo "- \`$OUT_DIR/massif/massif.ms_print.txt\`"
    fi
    echo
  } >> "$summary_path"
  if [[ -f "$OUT_DIR/massif/massif.ms_print.txt" ]]; then
    append_fenced_excerpt "$summary_path" "Massif Excerpt" "$OUT_DIR/massif/massif.ms_print.txt" 40
  fi
}

run_perf_stat() {
  local dir="$OUT_DIR/perf-stat"
  mkdir -p "$dir"
  log "running perf stat"
  "$PERF_BIN" stat -x ';' -o "$dir/perf.stat.csv" \
    -e "$PERF_STAT_EVENTS" -- \
    "$MAPPER_BIN" -I "$INDEX" -i "$READS_SUBSET" -o "$dir/mapper.sam" -k "$K" \
    > "$dir/mapper.stdout.log" 2> "$dir/mapper.stderr.log"
}

run_perf_record() {
  local dir="$OUT_DIR/perf-record"
  mkdir -p "$dir"
  log "running perf record"
  "$PERF_BIN" record -g -F "$PERF_RECORD_FREQ" -o "$dir/perf.data" -- \
    "$MAPPER_BIN" -I "$INDEX" -i "$READS_SUBSET" -o "$dir/mapper.sam" -k "$K" \
    > "$dir/mapper.stdout.log" 2> "$dir/mapper.stderr.log"
  log "rendering perf report"
  "$PERF_BIN" report --stdio -i "$dir/perf.data" > "$dir/perf.report.txt"
}

prepare_valgrind_mapper_bin() {
  local dir="$OUT_DIR/valgrind-bin"
  local backup_bin="$dir/mapper.native.backup"
  local portable_bin="$dir/mapper.portable"

  mkdir -p "$dir"
  if [[ "$VALGRIND_SAFE_BUILD" != "1" ]]; then
    VALGRIND_MAPPER_BIN="$MAPPER_BIN"
    return
  fi

  log "building portable mapper for valgrind"
  cp "$MAPPER_BIN" "$backup_bin"
  make -C "$ROOT_DIR" portable >/dev/null
  cp "$ROOT_DIR/mapper" "$portable_bin"
  cp "$backup_bin" "$ROOT_DIR/mapper"
  chmod +x "$portable_bin"
  VALGRIND_MAPPER_BIN="$portable_bin"
}

run_valgrind_command() {
  local tool_name="$1"
  shift
  local rc=0
  set +e
  env MAPPER_SPEED_MAX_SIMD="$VALGRIND_SIMD_CAP" "$@"
  rc=$?
  set -e
  if [[ "$rc" -eq 132 ]]; then
    die "$tool_name failed with illegal instruction. Try VALGRIND_SIMD_CAP=generic and VALGRIND_SAFE_BUILD=1 (defaults), or rebuild manually with 'make portable'."
  fi
  if [[ "$rc" -ne 0 ]]; then
    die "$tool_name failed with exit code $rc"
  fi
}

run_callgrind() {
  local dir="$OUT_DIR/callgrind"
  local -a callgrind_extra=()
  mkdir -p "$dir"
  read -r -a callgrind_extra <<< "$CALLGRIND_ARGS"
  log "running valgrind callgrind"
  run_valgrind_command "callgrind" \
    "$VALGRIND_BIN" --tool=callgrind --callgrind-out-file="$dir/callgrind.out" \
    "${callgrind_extra[@]}" \
    "$VALGRIND_MAPPER_BIN" -I "$INDEX" -i "$READS_SUBSET" -o "$dir/mapper.sam" -k "$K" \
    > "$dir/mapper.stdout.log" 2> "$dir/mapper.stderr.log"
  if command -v "$CALLGRIND_ANNOTATE_BIN" >/dev/null 2>&1; then
    log "rendering callgrind annotate"
    "$CALLGRIND_ANNOTATE_BIN" --auto=yes "$dir/callgrind.out" > "$dir/callgrind.annotate.txt"
  fi
}

run_massif() {
  local dir="$OUT_DIR/massif"
  local -a massif_extra=()
  mkdir -p "$dir"
  read -r -a massif_extra <<< "$MASSIF_ARGS"
  log "running valgrind massif"
  run_valgrind_command "massif" \
    "$VALGRIND_BIN" --tool=massif --massif-out-file="$dir/massif.out" \
    "${massif_extra[@]}" \
    "$VALGRIND_MAPPER_BIN" -I "$INDEX" -i "$READS_SUBSET" -o "$dir/mapper.sam" -k "$K" \
    > "$dir/mapper.stdout.log" 2> "$dir/mapper.stderr.log"
  if command -v "$MS_PRINT_BIN" >/dev/null 2>&1; then
    log "rendering ms_print"
    "$MS_PRINT_BIN" "$dir/massif.out" > "$dir/massif.ms_print.txt"
  fi
}

READS_SUBSET="$OUT_DIR/reads_subset.fastq"
extract_subset "$READS" "$PROFILE_READS" "$READS_SUBSET"
PROFILE_READ_COUNT="$(count_reads "$READS_SUBSET")"

SUMMARY_MD="$OUT_DIR/summary.md"
write_metadata "$OUT_DIR/run.meta"
append_markdown_header "$SUMMARY_MD"

{
  echo "Command template:"
  echo
  echo '```bash'
  echo "\"$MAPPER_BIN\" -I \"$INDEX\" -i \"$READS_SUBSET\" -o <tool-specific-output.sam> -k \"$K\""
  echo '```'
  echo
} >> "$SUMMARY_MD"

if mode_enabled perf-stat; then
  if command -v "$PERF_BIN" >/dev/null 2>&1; then
    run_perf_stat
    perf_stat_section "$OUT_DIR/perf-stat/perf.stat.csv" "$SUMMARY_MD"
  else
    append_tool_unavailable "$SUMMARY_MD" "Perf Stat" "Skipped because \`$PERF_BIN\` was not found in \`PATH\`."
  fi
fi

if mode_enabled perf-record; then
  if command -v "$PERF_BIN" >/dev/null 2>&1; then
    run_perf_record
    perf_record_section "$OUT_DIR/perf-record/perf.report.txt" "$SUMMARY_MD"
  else
    append_tool_unavailable "$SUMMARY_MD" "Perf Record" "Skipped because \`$PERF_BIN\` was not found in \`PATH\`."
  fi
fi

if mode_enabled callgrind; then
  if command -v "$VALGRIND_BIN" >/dev/null 2>&1; then
    if [[ -z "$VALGRIND_MAPPER_BIN" ]]; then
      prepare_valgrind_mapper_bin
    fi
    run_callgrind
    if [[ -f "$OUT_DIR/callgrind/callgrind.annotate.txt" ]]; then
      callgrind_section "$OUT_DIR/callgrind/callgrind.annotate.txt" "$SUMMARY_MD"
    else
      append_tool_unavailable "$SUMMARY_MD" "Callgrind" "Callgrind raw output was generated, but \`$CALLGRIND_ANNOTATE_BIN\` was not available to render an annotated report."
    fi
  else
    append_tool_unavailable "$SUMMARY_MD" "Callgrind" "Skipped because \`$VALGRIND_BIN\` was not found in \`PATH\`."
  fi
fi

if mode_enabled massif; then
  if command -v "$VALGRIND_BIN" >/dev/null 2>&1; then
    if [[ -z "$VALGRIND_MAPPER_BIN" ]]; then
      prepare_valgrind_mapper_bin
    fi
    run_massif
    massif_section "$OUT_DIR/massif/massif.out" "$SUMMARY_MD"
  else
    append_tool_unavailable "$SUMMARY_MD" "Massif" "Skipped because \`$VALGRIND_BIN\` was not found in \`PATH\`."
  fi
fi

log "summary written to $SUMMARY_MD"

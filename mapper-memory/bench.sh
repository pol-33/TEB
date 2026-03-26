#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$ROOT_DIR/.." && pwd)"

REF="${REF:-$REPO_DIR/data/genome.fa}"
READS="${READS:-$REPO_DIR/data/reads_1M.fastq}"
INDEX="${INDEX:-$ROOT_DIR/bench-results/genome.fmidx}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/bench-results}"
K="${K:-1}"
BENCH_READS="${BENCH_READS:-1000}"
CORRECTNESS_READS="${CORRECTNESS_READS:-250}"
BUILD_INDEX_IF_MISSING="${BUILD_INDEX_IF_MISSING:-0}"
TIME_STYLE="${TIME_STYLE:-auto}"

mkdir -p "$OUT_DIR"

log() {
  printf '[bench] %s\n' "$*"
}

warn() {
  printf '[bench][warn] %s\n' "$*" >&2
}

die() {
  printf '[bench][error] %s\n' "$*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: ./bench.sh

Environment overrides:
  REF=...                     Reference FASTA (default: ../data/genome.fa)
  READS=...                   FASTQ reads (default: ../data/reads_1M.fastq)
  INDEX=...                   FM-index output / input path
  OUT_DIR=...                 Benchmark output directory
  K=1                         Maximum edit distance (0..3)
  BENCH_READS=1000            Reads used for the timed mapper run; set to 0 for the full file
  CORRECTNESS_READS=250       Reads used for the optional correctness subset
  BUILD_INDEX_IF_MISSING=0    Set to 1 to build the FM-index automatically when INDEX is absent
  TIME_STYLE=auto             auto | full | portable
                              full tries '/usr/bin/time -l'; portable uses '/usr/bin/time -p'

Notes:
  - The contest judges mapping time and peak RSS excluding index build time.
  - This script records both index-build metrics and mapping metrics separately.
  - On a full GRCh38 reference, automatic index construction is a large one-time job.
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

[[ -f "$REF" ]] || die "reference FASTA not found: $REF"
[[ -f "$READS" ]] || die "reads FASTQ not found: $READS"
[[ "$K" =~ ^[0-3]$ ]] || die "K must be an integer in [0,3]"

count_reads() {
  awk 'END { printf "%d\n", NR / 4 }' "$1"
}

extract_subset() {
  local input_path="$1"
  local read_count="$2"
  local output_path="$3"
  if [[ "$read_count" -le 0 ]]; then
    cp "$input_path" "$output_path"
    return
  fi
  awk -v limit="$((read_count * 4))" 'NR <= limit { print }' "$input_path" > "$output_path"
}

has_bwa_index() {
  [[ -f "${REF}.amb" && -f "${REF}.ann" && -f "${REF}.bwt" && -f "${REF}.pac" && -f "${REF}.sa" ]]
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
    printf '%s\n' "full"
    return
  fi
  rm -f "$probe"
  printf '%s\n' "portable"
}

parse_peak_rss_from_stderr() {
  local stderr_path="$1"
  awk '
    /peak RSS:/ {
      line = $0;
      sub(/^.*peak RSS: /, "", line);
      sub(/ MiB.*$/, "", line);
      peak = int((line * 1024 * 1024) + 0.5);
    }
    /peak RSS / {
      line = $0;
      sub(/^.*peak RSS /, "", line);
      sub(/ MiB.*$/, "", line);
      peak = int((line * 1024 * 1024) + 0.5);
    }
    END {
      if (peak == "") peak = 0;
      print peak;
    }
  ' "$stderr_path"
}

run_with_time() {
  local label="$1"
  local stdout_path="$2"
  local stderr_path="$3"
  local time_path="$4"
  shift 4

  log "running $label"
  if [[ "$RESOLVED_TIME_STYLE" == "full" ]]; then
    /usr/bin/time -l -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
  else
    /usr/bin/time -p -o "$time_path" "$@" > "$stdout_path" 2> "$stderr_path"
  fi
  local rc=$?
  if [[ -f "$stderr_path" && -s "$stderr_path" ]]; then
    cat "$stderr_path" >&2
  fi
  return "$rc"
}

parse_time_report() {
  local time_path="$1"
  awk '
    /real/ && /user/ && /sys/            { real=$1; user=$3; sys=$5 }
    /^real[[:space:]]+/                  { real=$2 }
    /^user[[:space:]]+/                  { user=$2 }
    /^sys[[:space:]]+/                   { sys=$2 }
    /instructions retired/               { instr=$1 }
    /cycles elapsed/                     { cycles=$1 }
    /maximum resident set size/          { rss=$1 }
    /page faults/                        { pf=$1 }
    END {
      if (real == "") real = 0;
      if (user == "") user = 0;
      if (sys == "") sys = 0;
      if (instr == "") instr = 0;
      if (cycles == "") cycles = 0;
      if (rss == "") rss = 0;
      if (pf == "") pf = 0;
      printf "real=%s\nuser=%s\nsys=%s\ninstr=%s\ncycles=%s\nrss=%s\npf=%s\n", real, user, sys, instr, cycles, rss, pf;
    }
  ' "$time_path"
}

write_chrom_sizes() {
  local ref_path="$1"
  local output_path="$2"
  awk '
    /^>/ {
      if (name != "") {
        printf "%s\t%d\n", name, len;
      }
      name = substr($0, 2);
      sub(/[[:space:]].*$/, "", name);
      len = 0;
      next;
    }
    {
      gsub(/[[:space:]]/, "", $0);
      len += length($0);
    }
    END {
      if (name != "") {
        printf "%s\t%d\n", name, len;
      }
    }
  ' "$ref_path" > "$output_path"
}

write_failure_summary() {
  local phase="$1"
  local reason="$2"
  local summary_path="$OUT_DIR/summary.txt"

  {
    echo "Benchmark Summary"
    echo "================="
    echo "reference: $REF"
    echo "reads: $READS"
    echo "index: $INDEX"
    echo "k: $K"
    echo "timing_backend: ${RESOLVED_TIME_STYLE:-unknown}"
    echo
    echo "status: FAILED"
    echo "phase: $phase"
    echo "reason: $reason"
    echo
    if [[ -f "$OUT_DIR/${phase}.metrics" ]]; then
      echo "Captured Metrics"
      echo "----------------"
      cat "$OUT_DIR/${phase}.metrics"
    fi
  } > "$summary_path"

  cat "$summary_path"
}

validate_mapper_output() {
  local output_path="$1"
  local chrom_sizes_path="$2"
  local expected_reads="$3"
  local report_path="$4"

  awk -v expected_reads="$expected_reads" '
    function cigar_query_len(cigar,    rest, n, op, len) {
      rest = cigar;
      len = 0;
      while (match(rest, /^[0-9]+[MID]/)) {
        n = substr(rest, RSTART, RLENGTH - 1) + 0;
        op = substr(rest, RSTART + RLENGTH - 1, 1);
        if (op == "M" || op == "I") {
          len += n;
        }
        rest = substr(rest, RSTART + RLENGTH);
      }
      return (rest == "" ? len : -1);
    }
    function cigar_ref_len(cigar,    rest, n, op, len) {
      rest = cigar;
      len = 0;
      while (match(rest, /^[0-9]+[MID]/)) {
        n = substr(rest, RSTART, RLENGTH - 1) + 0;
        op = substr(rest, RSTART + RLENGTH - 1, 1);
        if (op == "M" || op == "D") {
          len += n;
        }
        rest = substr(rest, RSTART + RLENGTH);
      }
      return (rest == "" ? len : -1);
    }
    function fail(msg) {
      print msg > "/dev/stderr";
      exit 1;
    }
    BEGIN {
      FS = "[ \t]+";
    }
    FNR == NR {
      chrom_len[$1] = $2 + 0;
      next;
    }
    {
      total++;
      if (NF < 6 || NF > 7) {
        fail("invalid field count at output line " total);
      }
      read_name = $1;
      chrom = $2;
      pos = $3 + 0;
      cigar = $4;
      seq = $5;
      qual = $6;
      alt = (NF == 7 ? $7 : "");

      if (length(seq) != length(qual)) {
        fail("sequence/quality length mismatch for " read_name);
      }

      if (chrom == "*") {
        if ($3 != "0" || cigar != "*") {
          fail("unmapped record has invalid coordinates for " read_name);
        }
        unmapped++;
      } else {
        if (!(chrom in chrom_len)) {
          fail("unknown chromosome " chrom " for " read_name);
        }
        if (pos < 1 || pos > chrom_len[chrom]) {
          fail("position out of bounds for " read_name);
        }
        if (cigar !~ /^([0-9]+[MID])+$/) {
          fail("invalid primary cigar for " read_name ": " cigar);
        }
        qlen = cigar_query_len(cigar);
        rlen = cigar_ref_len(cigar);
        if (qlen != length(seq)) {
          fail("primary cigar/query length mismatch for " read_name);
        }
        if (rlen < 1 || pos + rlen - 1 > chrom_len[chrom]) {
          fail("primary alignment crosses chromosome end for " read_name);
        }
        mapped++;
      }

      if (alt != "") {
        if (alt !~ /^ALT:[^,]+,[0-9]+,([0-9]+[MID])+$/) {
          fail("invalid ALT field for " read_name ": " alt);
        }
        sub(/^ALT:/, "", alt);
        split(alt, alt_parts, ",");
        alt_chrom = alt_parts[1];
        alt_pos = alt_parts[2] + 0;
        alt_cigar = alt_parts[3];
        if (!(alt_chrom in chrom_len)) {
          fail("unknown ALT chromosome for " read_name);
        }
        alt_qlen = cigar_query_len(alt_cigar);
        alt_rlen = cigar_ref_len(alt_cigar);
        if (alt_qlen != length(seq)) {
          fail("ALT cigar/query length mismatch for " read_name);
        }
        if (alt_pos < 1 || alt_pos + alt_rlen - 1 > chrom_len[alt_chrom]) {
          fail("ALT alignment crosses chromosome end for " read_name);
        }
        alt_count++;
      }
    }
    END {
      if (total != expected_reads) {
        fail("output line count mismatch: expected " expected_reads ", got " total);
      }
      printf "reads=%d\nmapped=%d\nunmapped=%d\nalt=%d\n", total, mapped, unmapped, alt_count;
    }
  ' "$chrom_sizes_path" "$output_path" > "$report_path"
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
    function flag_set(flag, mask) {
      return int(flag / mask) % 2;
    }
    function has_soft_or_hard_clip(cigar) {
      return cigar ~ /[SH]/;
    }
    /^@/ { next }
    {
      flag = $2 + 0;
      if (flag_set(flag, 256) || flag_set(flag, 2048)) {
        next;
      }
      if ($3 == "*" || flag_set(flag, 4)) {
        print $1 "\t*\t0\t*";
        next;
      }
      nm = -1;
      for (i = 12; i <= NF; ++i) {
        if ($i ~ /^NM:i:/) {
          nm = substr($i, 6) + 0;
        }
      }
      if (nm < 0 || nm > max_nm || has_soft_or_hard_clip($6)) {
        print $1 "\t*\t0\t*";
        next;
      }
      print $1 "\t" $3 "\t" $4 "\t" $6;
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
      mapper_is_mapped = (name in mapper_chr) && mapper_chr[name] != "*";
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
          if (delta <= 50) within50++;
          if (mapper_cigar[name] == bwa_cigar) exact_cigar++;
        }
      } else if (mapper_is_mapped && !bwa_is_mapped) {
        mapper_only++;
      } else if (!mapper_is_mapped && bwa_is_mapped) {
        bwa_only++;
      }
    }
    END {
      printf "subset_reads=%d\n", total;
      printf "mapper_mapped=%d\n", mapper_mapped;
      printf "bwa_mapped=%d\n", bwa_mapped;
      printf "shared_mapped=%d\n", shared;
      printf "chrom_match=%d\n", chrom_match;
      printf "within_10bp=%d\n", within10;
      printf "within_50bp=%d\n", within50;
      printf "exact_pos=%d\n", exact_pos;
      printf "exact_cigar=%d\n", exact_cigar;
      printf "mapper_only=%d\n", mapper_only;
      printf "bwa_only=%d\n", bwa_only;
    }
  ' "$mapper_tsv" "$bwa_tsv" > "$report_path"
}

log "building mapper-memory binaries"
make -C "$ROOT_DIR"
RESOLVED_TIME_STYLE="$(detect_time_style)"
log "timing backend: $RESOLVED_TIME_STYLE"

CHROM_SIZES="$OUT_DIR/chrom.sizes.tsv"
if [[ -s "$CHROM_SIZES" && "$CHROM_SIZES" -nt "$REF" ]]; then
  log "using cached chromosome sizes: $CHROM_SIZES"
else
  log "scanning reference to compute chromosome sizes"
  write_chrom_sizes "$REF" "$CHROM_SIZES"
fi

BENCH_FASTQ="$OUT_DIR/bench_reads.fastq"
log "preparing benchmark read set"
extract_subset "$READS" "$BENCH_READS" "$BENCH_FASTQ"
BENCH_READ_COUNT="$(count_reads "$BENCH_FASTQ")"

if [[ -f "$INDEX" ]]; then
  log "using existing FM-index: $INDEX"
else
  if [[ "$BUILD_INDEX_IF_MISSING" != "1" ]]; then
    warn "FM-index not found at $INDEX"
    warn "Set BUILD_INDEX_IF_MISSING=1 to time index construction automatically."
  else
    log "full FM-index build requested; on GRCh38 this can take a long time before mapping starts"
	    INDEX_STDOUT="$OUT_DIR/index.stdout.log"
	    INDEX_STDERR="$OUT_DIR/index.stderr.log"
	    INDEX_TIME="$OUT_DIR/index.time.log"
	    if run_with_time "indexer" "$INDEX_STDOUT" "$INDEX_STDERR" "$INDEX_TIME" "$ROOT_DIR/indexer" -R "$REF" -I "$INDEX"; then
	      parse_time_report "$INDEX_TIME" > "$OUT_DIR/index.metrics"
        if [[ "$RESOLVED_TIME_STYLE" != "full" ]]; then
          index_peak_rss="$(parse_peak_rss_from_stderr "$INDEX_STDERR")"
          awk -v rss_override="$index_peak_rss" '
            BEGIN { updated = 0 }
            /^rss=/ {
              print "rss=" rss_override;
              updated = 1;
              next;
            }
            { print }
            END {
              if (!updated) {
                print "rss=" rss_override;
              }
            }
          ' "$OUT_DIR/index.metrics" > "$OUT_DIR/index.metrics.tmp"
          mv "$OUT_DIR/index.metrics.tmp" "$OUT_DIR/index.metrics"
        fi
	    else
	      if [[ -f "$INDEX_TIME" ]]; then
	        parse_time_report "$INDEX_TIME" > "$OUT_DIR/index.metrics"
          if [[ "$RESOLVED_TIME_STYLE" != "full" ]]; then
            index_peak_rss="$(parse_peak_rss_from_stderr "$INDEX_STDERR")"
            awk -v rss_override="$index_peak_rss" '
              BEGIN { updated = 0 }
              /^rss=/ {
                print "rss=" rss_override;
                updated = 1;
                next;
              }
              { print }
              END {
                if (!updated) {
                  print "rss=" rss_override;
                }
              }
            ' "$OUT_DIR/index.metrics" > "$OUT_DIR/index.metrics.tmp"
            mv "$OUT_DIR/index.metrics.tmp" "$OUT_DIR/index.metrics"
          fi
	      fi
	      write_failure_summary "index" "FM-index construction failed before mapping. See $INDEX_STDERR"
	      exit 1
	    fi
  fi
fi

if [[ ! -f "$INDEX" ]]; then
  warn "skipping mapper benchmark because no FM-index is available."
  exit 0
fi

MAP_OUTPUT="$OUT_DIR/mapper.out.sam"
MAP_TIME="$OUT_DIR/mapper.time.log"
MAP_STDERR="$OUT_DIR/mapper.stderr.log"

log "benchmarking mapper on $BENCH_READ_COUNT read(s) with k=$K"
run_with_time "mapper" /dev/null "$MAP_STDERR" "$MAP_TIME" "$ROOT_DIR/mapper" -I "$INDEX" -i "$BENCH_FASTQ" -o "$MAP_OUTPUT" -k "$K" || {
  cat "$MAP_STDERR" >&2
  die "mapper run failed"
}

parse_time_report "$MAP_TIME" > "$OUT_DIR/mapper.metrics"
if [[ "$RESOLVED_TIME_STYLE" != "full" ]]; then
  mapper_peak_rss="$(parse_peak_rss_from_stderr "$MAP_STDERR")"
  awk -v rss_override="$mapper_peak_rss" '
    BEGIN { updated = 0 }
    /^rss=/ {
      print "rss=" rss_override;
      updated = 1;
      next;
    }
    { print }
    END {
      if (!updated) {
        print "rss=" rss_override;
      }
    }
  ' "$OUT_DIR/mapper.metrics" > "$OUT_DIR/mapper.metrics.tmp"
  mv "$OUT_DIR/mapper.metrics.tmp" "$OUT_DIR/mapper.metrics"
fi
validate_mapper_output "$MAP_OUTPUT" "$CHROM_SIZES" "$BENCH_READ_COUNT" "$OUT_DIR/mapper.validation"

source "$OUT_DIR/mapper.metrics"
source "$OUT_DIR/mapper.validation"

{
  echo "Benchmark Summary"
  echo "================="
  echo "reference: $REF"
  echo "reads: $READS"
  echo "bench_reads: $BENCH_READ_COUNT"
  echo "index: $INDEX"
  echo "k: $K"
  echo "timing_backend: $RESOLVED_TIME_STYLE"
  echo
  echo "Mapping Metrics"
  echo "---------------"
  echo "real_seconds: $real"
  echo "user_seconds: $user"
  echo "sys_seconds: $sys"
  echo "peak_rss_bytes: $rss"
  echo "page_faults: $pf"
  echo "instructions: $instr"
  echo "cycles: $cycles"
  awk -v reads="$BENCH_READ_COUNT" -v real="$real" 'BEGIN { if (real > 0) printf "reads_per_second: %.2f\n", reads / real; else print "reads_per_second: 0" }'
  echo
  echo "Output Validation"
  echo "-----------------"
  echo "mapped: $mapped"
  echo "unmapped: $unmapped"
  echo "with_alt: $alt"
} > "$OUT_DIR/summary.txt"

cat "$OUT_DIR/summary.txt"

if [[ "$RUN_BWA_CHECK" != "1" ]]; then
  log "skipping bwa correctness comparison by request"
  exit 0
fi

if ! command -v bwa >/dev/null 2>&1; then
  warn "bwa is not installed; skipping professional-baseline correctness check"
  exit 0
fi

if ! has_bwa_index; then
  if [[ "$BUILD_BWA_INDEX_IF_MISSING" != "1" ]]; then
    warn "bwa index files are missing for $REF"
    warn "Set BUILD_BWA_INDEX_IF_MISSING=1 to build them automatically."
    exit 0
  fi
  BWA_INDEX_STDOUT="$OUT_DIR/bwa-index.stdout.log"
  BWA_INDEX_STDERR="$OUT_DIR/bwa-index.stderr.log"
  BWA_INDEX_TIME="$OUT_DIR/bwa-index.time.log"
  run_with_time "bwa index" "$BWA_INDEX_STDOUT" "$BWA_INDEX_STDERR" "$BWA_INDEX_TIME" bwa index "$REF"
  parse_time_report "$BWA_INDEX_TIME" > "$OUT_DIR/bwa-index.metrics"
fi

CORRECTNESS_FASTQ="$OUT_DIR/correctness_reads.fastq"
extract_subset "$READS" "$CORRECTNESS_READS" "$CORRECTNESS_FASTQ"
CORRECTNESS_COUNT="$(count_reads "$CORRECTNESS_FASTQ")"

MAPPER_CORRECTNESS_OUTPUT="$OUT_DIR/correctness.mapper.sam"
log "running mapper correctness subset ($CORRECTNESS_COUNT reads)"
"$ROOT_DIR/mapper" -I "$INDEX" -i "$CORRECTNESS_FASTQ" -o "$MAPPER_CORRECTNESS_OUTPUT" -k "$K" > /dev/null
validate_mapper_output "$MAPPER_CORRECTNESS_OUTPUT" "$CHROM_SIZES" "$CORRECTNESS_COUNT" "$OUT_DIR/correctness.validation"

BWA_SAM="$OUT_DIR/correctness.bwa.sam"
BWA_STDERR="$OUT_DIR/correctness.bwa.stderr.log"
BWA_TIME="$OUT_DIR/correctness.bwa.time.log"
run_with_time "bwa mem" "$BWA_SAM" "$BWA_STDERR" "$BWA_TIME" bwa mem "$REF" "$CORRECTNESS_FASTQ"
parse_time_report "$BWA_TIME" > "$OUT_DIR/correctness.bwa.metrics"

MAPPER_TSV="$OUT_DIR/correctness.mapper.tsv"
BWA_TSV="$OUT_DIR/correctness.bwa.tsv"
convert_mapper_output "$MAPPER_CORRECTNESS_OUTPUT" "$MAPPER_TSV"
convert_bwa_output "$BWA_SAM" "$BWA_TSV" "$K"
compare_against_bwa "$MAPPER_TSV" "$BWA_TSV" "$OUT_DIR/correctness.compare"

{
  echo
  echo "Correctness Against bwa mem"
  echo "---------------------------"
  cat "$OUT_DIR/correctness.compare"
} >> "$OUT_DIR/summary.txt"

cat "$OUT_DIR/correctness.compare"

#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$ROOT_DIR/.." && pwd)"

REF="${REF:-$REPO_DIR/data/genome.fa}"
READS="${READS:-$REPO_DIR/data/reads_1M.fastq}"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/bench-small}"
SAMPLE_READS="${SAMPLE_READS:-200000}"
SELECTED_READS="${SELECTED_READS:-128}"
WINDOW_BEFORE="${WINDOW_BEFORE:-250}"
WINDOW_AFTER="${WINDOW_AFTER:-350}"
MAPQ_MIN="${MAPQ_MIN:-60}"
MIN_GAP="${MIN_GAP:-5000}"
BWA_THREADS="${BWA_THREADS:-1}"

READ_SAMPLE="$OUT_DIR/reads_sample.fastq"
ALIGNMENT_SAM="$OUT_DIR/fullgenome.sample.sam"
LOCI_TSV="$OUT_DIR/selected_loci.tsv"
REFERENCE_OUT="$OUT_DIR/reference.fa"
READ_IDS="$OUT_DIR/selected_read_ids.txt"
READS_OUT="$OUT_DIR/reads.fastq"
SUMMARY_OUT="$OUT_DIR/summary.txt"

mkdir -p "$OUT_DIR"

extract_region() {
  local chr="$1"
  local start="$2"
  local end="$3"
  local name="$4"
  local fai_line
  fai_line="$(awk -F'\t' -v chr="$chr" '$1==chr {print $0}' "${REF}.fai")"
  [[ -n "$fai_line" ]] || { echo "missing $chr in ${REF}.fai" >&2; exit 1; }

  local chrom_len file_offset line_bases line_bytes
  chrom_len="$(printf '%s\n' "$fai_line" | awk -F'\t' '{print $2}')"
  file_offset="$(printf '%s\n' "$fai_line" | awk -F'\t' '{print $3}')"
  line_bases="$(printf '%s\n' "$fai_line" | awk -F'\t' '{print $4}')"
  line_bytes="$(printf '%s\n' "$fai_line" | awk -F'\t' '{print $5}')"

  if (( end > chrom_len )); then
    end="$chrom_len"
  fi

  local start0=$((start - 1))
  local end0=$((end - 1))
  local byte_start=$((file_offset + (start0 / line_bases) * line_bytes + (start0 % line_bases)))
  local byte_end=$((file_offset + (end0 / line_bases) * line_bytes + (end0 % line_bases)))
  local byte_count=$((byte_end - byte_start + 1))

  {
    printf '>%s\n' "$name"
    (
      set +o pipefail
      tail -c +"$((byte_start + 1))" "$REF" | head -c "$byte_count" | tr -d '\n' | fold -w 60
    )
    printf '\n'
  } >> "$REFERENCE_OUT"
}

echo "[prepare-small] taking first $SAMPLE_READS reads from $READS"
awk -v limit="$((SAMPLE_READS * 4))" 'NR <= limit { print }' "$READS" > "$READ_SAMPLE"

echo "[prepare-small] aligning sample reads to full genome with bwa mem"
bwa mem -t "$BWA_THREADS" "$REF" "$READ_SAMPLE" > "$ALIGNMENT_SAM"

echo "[prepare-small] selecting up to $SELECTED_READS high-confidence loci"
awk -v fai_path="${REF}.fai" \
    -v selected_reads="$SELECTED_READS" \
    -v window_before="$WINDOW_BEFORE" \
    -v window_after="$WINDOW_AFTER" \
    -v mapq_min="$MAPQ_MIN" \
    -v min_gap="$MIN_GAP" '
  BEGIN {
    while ((getline line < fai_path) > 0) {
      split(line, f, "\t");
      chrom_len[f[1]] = f[2] + 0;
    }
  }
  /^@/ { next }
  count >= selected_reads { next }
  $3 == "*" { next }
  index($3, "_") > 0 { next }
  $5 + 0 < mapq_min { next }
  {
    bucket = $3 ":" int(($4 + 0) / min_gap);
    if (bucket in seen_bucket) {
      next;
    }
    start = $4 - window_before;
    if (start < 1) {
      start = 1;
    }
    end = $4 + window_after;
    if (end > chrom_len[$3]) {
      end = chrom_len[$3];
    }
    print $1 "\t" $3 "\t" start "\t" end;
    seen_bucket[bucket] = 1;
    ++count;
  }
' "$ALIGNMENT_SAM" > "$LOCI_TSV"

LOCUS_COUNT="$(wc -l < "$LOCI_TSV" | tr -d ' ')"
if [[ "$LOCUS_COUNT" -eq 0 ]]; then
  echo "no loci selected" >&2
  exit 1
fi

echo "[prepare-small] extracting reference windows to $REFERENCE_OUT"
: > "$REFERENCE_OUT"
awk -F'\t' '{print $1}' "$LOCI_TSV" > "$READ_IDS"
while IFS=$'\t' read -r qname chr start end; do
  extract_region "$chr" "$start" "$end" "${chr}_${start}_${end}"
done < "$LOCI_TSV"

echo "[prepare-small] extracting selected reads to $READS_OUT"
awk -v ids_path="$READ_IDS" '
  BEGIN {
    while ((getline line < ids_path) > 0) {
      wanted[line] = 1;
    }
  }
  NR % 4 == 1 {
    header = substr($0, 2);
    sub(/[ \t].*$/, "", header);
    keep = (header in wanted);
  }
  keep { print }
' "$READ_SAMPLE" > "$READS_OUT"

SELECTED_COUNT="$(awk 'END { printf "%d\n", NR / 4 }' "$READS_OUT")"
TOTAL_BP="$(awk -F'\t' '{sum += ($4 - $3 + 1)} END {print sum + 0}' "$LOCI_TSV")"

{
  echo "Small Benchmark Dataset"
  echo "======================="
  echo "reference: $REFERENCE_OUT"
  echo "reads: $READS_OUT"
  echo "sample_reads_scanned: $SAMPLE_READS"
  echo "selected_reads: $SELECTED_COUNT"
  echo "loci: $LOCUS_COUNT"
  echo "reference_total_bp: $TOTAL_BP"
  echo "mapq_min: $MAPQ_MIN"
  echo "window_before: $WINDOW_BEFORE"
  echo "window_after: $WINDOW_AFTER"
  echo "min_gap: $MIN_GAP"
  echo
  echo "Selected Loci"
  echo "-------------"
  cat "$LOCI_TSV"
} | tee "$SUMMARY_OUT"

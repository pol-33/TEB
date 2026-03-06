#!/usr/bin/env bash
# =============================================================================
#  bench.sh  —  TEB parser benchmark across git commits
#
#  This machine is macOS / Apple Silicon — Linux's "perf" is not available.
#  We use  /usr/bin/time -l  instead, which exposes the same hardware counters
#  through the kernel's performance monitoring unit (PMU):
#
#    Metric reported by /usr/bin/time -l   │  Linux perf stat equivalent
#    ──────────────────────────────────────┼──────────────────────────────────
#    real (wall clock)                     │  elapsed time
#    user + sys (CPU time)                 │  task-clock
#    instructions retired                  │  perf -e instructions
#    cycles elapsed                        │  perf -e cycles
#    maximum resident set size             │  peak RSS (/usr/bin/time on Linux)
#    page faults                           │  perf -e page-faults
#    IPC = instructions / cycles           │  perf stat (derived)
#
#  Runs:
#    FASTA  — 4 commits, FASTA_RUNS timed runs + 1 warm-up each
#    FASTQ  — 2 commits, FASTQ_RUNS timed runs + 1 warm-up each
#             ⚠ Pre-mmap FASTQ takes ~247 s / run — set FASTQ_RUNS=1 unless patient
#
#  Usage:  ./bench.sh
# =============================================================================
set -uo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA_IN="$REPO/datasets/chr1.fna"
FASTQ_IN="$REPO/datasets/SRR22320000_1.fastq"
FASTA_OUT="/tmp/teb_bench_out.fna"
FASTQ_OUT="/tmp/teb_bench_out.fastq"

# ── Tuning ────────────────────────────────────────────────────────────────────
FASTA_RUNS=3   # timed runs (a warm-up is always prepended and discarded)
FASTQ_RUNS=1   # keep at 1 unless you have ~8 min free per extra run

# ── Commits ───────────────────────────────────────────────────────────────────
# Format: "short_hash:Label shown in table"
FASTA_COMMITS=(
    "0126f47:Make output file optional; add README:long"
    "c30cffa:Updating Kmer table + exercises"
    "8441748:Merge fixes"
    "58cca05:Optimize parsers"
    "b3e0979:mmap + per-seq flag (HEAD)"
)
FASTQ_COMMITS=(
    "58cca05:Optimize parsers"
    "b3e0979:mmap + per-seq flag (HEAD)"
)

# ── Colours (disabled when stdout is not a terminal) ─────────────────────────
if [[ -t 1 ]]; then
    BOLD=$'\e[1m' DIM=$'\e[2m' CYN=$'\e[36m'
    GRN=$'\e[32m' YEL=$'\e[33m' RED=$'\e[31m' RST=$'\e[0m'
else
    BOLD='' DIM='' CYN='' GRN='' YEL='' RED='' RST=''
fi

# ─────────────────────────────────────────────────────────────────────────────
# single_run <out_path> <out_flag> <binary> [extra_args...]
#
#   Executes:  /usr/bin/time -l <binary> [extra_args...] <out_flag> <out_path>
#   Discards program stdout (stats lines) with >/dev/null.
#   Captures /usr/bin/time stderr to a temp file, then parses it with awk.
#   Prints one line to stdout:
#     real_s  user_s  sys_s  instructions  cycles  rss_bytes  page_faults
# ─────────────────────────────────────────────────────────────────────────────
single_run() {
    local out_path="$1" out_flag="$2"; shift 2
    local tmpf; tmpf=$(mktemp)

    # Redirect program stdout to /dev/null (parser stat printout is O(1) but
    # adds noise); capture /usr/bin/time report from stderr.
    { /usr/bin/time -l "$@" "$out_flag" "$out_path" > /dev/null ; } 2>"$tmpf"

    awk '
        # First line of /usr/bin/time -l: "  real_s real   user_s user   sys_s sys"
        /real/ && /user/ && /sys/ { real=$1; user=$3; sys=$5 }
        /instructions retired/    { instr=$1  }
        /cycles elapsed/          { cycles=$1 }
        /maximum resident set/    { rss=$1    }
        /page faults/             { pf=$1     }
        END { print real, user, sys, instr, cycles, rss, pf }
    ' "$tmpf"

    rm -f "$tmpf"
}

# ─────────────────────────────────────────────────────────────────────────────
# benchmark_suite <format> <nruns> <out_path> <commits_array_name>
# ─────────────────────────────────────────────────────────────────────────────
benchmark_suite() {
    local fmt="$1" nruns="$2" out_path="$3"
    # Copy the named array into a local one (eval needed for bash 3.2 compat;
    # macOS ships bash 3.2 which lacks nameref  local -n)
    local -a _commits
    eval "_commits=(\"\${${4}[@]}\")"

    local fmt_upper; fmt_upper=$(echo "$fmt" | tr '[:lower:]' '[:upper:]')
    local inp; [[ $fmt == fasta ]] && inp="$FASTA_IN" || inp="$FASTQ_IN"
    local sep; sep=$(printf '─%.0s' {1..104})

    printf "\n%s\n" "$sep"
    printf "  %s%s%s PARSER — %d timed run(s) + 1 warm-up   │   input: %s\n" \
        "$BOLD$CYN" "$fmt_upper" "$RST" "$nruns" "$inp"
    printf "%s\n" "$sep"
    printf "%-32s  %8s %8s %8s  %12s %12s  %8s  %5s  %10s  %s\n" \
        "Commit" "Min(s)" "Avg(s)" "CPU(s)" "Instr(M)" "Cycles(M)" \
        "RSS(MB)" "IPC" "PgFaults" "Speedup"
    printf "%s\n" "$sep"

    # ── Save current git state ────────────────────────────────────────────────
    local saved_ref stashed=0
    saved_ref=$(git -C "$REPO" symbolic-ref --short HEAD 2>/dev/null \
                || git -C "$REPO" rev-parse HEAD)

    # Stash any dirty working tree so checkout never complains
    if ! git -C "$REPO" diff --quiet HEAD 2>/dev/null; then
        git -C "$REPO" stash push -m "bench.sh auto-stash" --quiet
        stashed=1
    fi

    local baseline_min=""   # min wall time of the first commit (for speedup column)
    local baseline_label="" # label of the first commit

    for entry in "${_commits[@]}"; do
        local hash="${entry%%:*}"
        local rest="${entry#*:}"      # "label" or "label:long"
        local label="${rest%%:*}"     # always the label
        local flags_style=""
        [[ "$rest" == *:* ]] && flags_style="${rest#*:}"

        # Select short or long CLI flags based on per-entry annotation
        local f_in="-i" f_fmt="-f" f_out="-o"
        if [[ "${flags_style:-}" == "long" ]]; then
            f_in="--input"; f_fmt="--format"; f_out="--output"
        fi

        # ── Checkout & build ─────────────────────────────────────────────────
        printf "  %s● building  %s  →  %s%s\n" "$YEL" "$hash" "$label" "$RST" >&2
        git -C "$REPO" checkout --quiet "$hash"
        rm -f "$REPO/teb" "$REPO/teb.exe"

        local build_log; build_log=$(mktemp)
        if ! make -C "$REPO" re > "$build_log" 2>&1; then
            printf "  %s[BUILD FAILED for %s — skipping. See %s]%s\n" \
                "$RED" "$hash" "$build_log" "$RST" >&2
            continue
        fi
        rm -f "$build_log"

        local binary=""
        if   [[ -x "$REPO/teb.exe" ]]; then binary="$REPO/teb.exe"
        elif [[ -x "$REPO/teb"     ]]; then binary="$REPO/teb"
        else
            printf "  %s[NO BINARY found after build of %s — skipping]%s\n" \
                "$RED" "$hash" "$RST" >&2
            continue
        fi

        # ── Warm-up run (page-cache priming, not measured) ───────────────────
        printf "  %s  warm-up ...%s\r" "$DIM" "$RST" >&2
        "$binary" "$f_in" "$inp" "$f_fmt" "$fmt" "$f_out" "$out_path" > /dev/null 2>&1 || true

        # ── Timed runs ───────────────────────────────────────────────────────
        local sum_real=0 sum_cpu=0 sum_instr=0 sum_cycles=0 sum_pf=0
        local min_real=9999999 last_rss=0

        for ((r = 1; r <= nruns; r++)); do
            printf "  %s  run %d / %d  ...%s\r" "$DIM" "$r" "$nruns" "$RST" >&2

            local real user sys instr cycles rss pf
            read -r real user sys instr cycles rss pf < <(
                single_run "$out_path" "$f_out" "$binary" "$f_in" "$inp" "$f_fmt" "$fmt"
            ) || true

            local cpu; cpu=$(awk "BEGIN{printf \"%.6f\", $user + $sys}")
            min_real=$(awk "BEGIN{ print ($real < $min_real) ? $real : $min_real }")
            sum_real=$(awk   "BEGIN{printf \"%.6f\", $sum_real  + $real  }")
            sum_cpu=$(awk    "BEGIN{printf \"%.6f\", $sum_cpu   + $cpu   }")
            sum_instr=$(awk  "BEGIN{printf \"%.0f\", $sum_instr + $instr }")
            sum_cycles=$(awk "BEGIN{printf \"%.0f\", $sum_cycles+ $cycles}")
            sum_pf=$(awk     "BEGIN{printf \"%.0f\", $sum_pf    + $pf    }")
            last_rss=$rss
        done

        printf "%-60s\r" "" >&2  # erase progress line

        # ── Compute aggregates ───────────────────────────────────────────────
        local avg_real avg_cpu m_instr m_cycles rss_mb ipc avg_pf speedup
        avg_real=$(awk  "BEGIN{printf \"%.3f\", $sum_real  / $nruns      }")
        avg_cpu=$(awk   "BEGIN{printf \"%.3f\", $sum_cpu   / $nruns      }")
        m_instr=$(awk   "BEGIN{printf \"%.1f\", $sum_instr / $nruns / 1e6}")
        m_cycles=$(awk  "BEGIN{printf \"%.1f\", $sum_cycles/ $nruns / 1e6}")
        rss_mb=$(awk    "BEGIN{printf \"%.1f\", $last_rss  / 1048576     }")
        ipc=$(awk       "BEGIN{c=$sum_cycles; printf \"%.2f\", (c>0) ? $sum_instr/c : 0}")
        avg_pf=$(awk    "BEGIN{printf \"%.0f\", $sum_pf    / $nruns      }")

        # ── Speedup vs baseline ──────────────────────────────────────────────
        if [[ -z "$baseline_min" ]]; then
            baseline_min="$min_real"
            baseline_label="$label"
            speedup="${DIM}(baseline)${RST}"
        else
            local sx; sx=$(awk "BEGIN{printf \"%.1f\", $baseline_min / $min_real}")
            speedup="${GRN}${BOLD}${sx}×${RST}"
        fi

        printf "%-32s  %8.3f %8.3f %8.3f  %12s %12s  %8s  %5s  %10s  %s\n" \
            "${label:0:32}" "$min_real" "$avg_real" "$avg_cpu" \
            "$m_instr" "$m_cycles" "$rss_mb" "$ipc" "$avg_pf" "$speedup"

        rm -f "$out_path"
    done

    printf "%s\n" "$sep"
    printf "  %sMin = fastest of %d run(s).  Speedup relative to \"%s\".%s\n" \
        "$DIM" "$nruns" "$baseline_label" "$RST"
    printf "  %sCPU = user+sys.  IPC = Instructions / Cycles (higher → better pipeline use).%s\n" \
        "$DIM" "$RST"

    # ── Restore git state ─────────────────────────────────────────────────────
    printf "\n  %s● restoring branch: %s%s\n" "$YEL" "$saved_ref" "$RST" >&2
    git -C "$REPO" checkout --force --quiet "$saved_ref"
    rm -f "$REPO/teb" "$REPO/teb.exe"
    make -C "$REPO" re > /dev/null 2>&1

    if [[ $stashed -eq 1 ]]; then
        git -C "$REPO" stash pop --quiet 2>/dev/null || true
    fi
}

# ── Main ──────────────────────────────────────────────────────────────────────
printf "\n%s══════════════════════════════════════════════════════════════════%s\n" \
    "$BOLD" "$RST"
printf "%s  TEB Parser Benchmark — %s%s\n" \
    "$BOLD" "$(date '+%Y-%m-%d %H:%M:%S')" "$RST"
printf "%s  Platform : %s %s (%s)%s\n" \
    "$DIM" "$(sw_vers -productName 2>/dev/null || uname -s)" \
    "$(sw_vers -productVersion 2>/dev/null || uname -r)" "$(uname -m)" "$RST"
printf "%s  Compiler : $(g++ --version | head -1)%s\n" "$DIM" "$RST"
printf "%s  NOTE     : macOS has no 'perf'. Using /usr/bin/time -l%s\n" \
    "$DIM" "$RST"
printf "%s            (instructions retired & cycles elapsed come from the PMU)%s\n" \
    "$DIM" "$RST"
printf "%s══════════════════════════════════════════════════════════════════%s\n\n" \
    "$BOLD" "$RST"

printf "%s  FASTA: 5 commits × %d run(s) + 1 warm-up%s\n" "$DIM" "$FASTA_RUNS" "$RST"
printf "%s  FASTQ: 2 commits × %d run(s) + 1 warm-up  (pre-mmap run ≈ 247 s)%s\n\n" \
    "$DIM" "$FASTQ_RUNS" "$RST"

benchmark_suite "fasta" "$FASTA_RUNS" "$FASTA_OUT" FASTA_COMMITS
benchmark_suite "fastq" "$FASTQ_RUNS" "$FASTQ_OUT" FASTQ_COMMITS

printf "\n%sBenchmark complete.%s\n\n" "$GRN$BOLD" "$RST"

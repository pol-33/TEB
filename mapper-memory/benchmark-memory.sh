#!/bin/bash
# Benchmark script to measure memory usage across optimization levels
# Usage: ./benchmark-memory.sh <genome.fa> <reads.fastq> [num_reads]

set -e

GENOME="$1"
READS="$2"
NUM_READS="${3:-100000}"  # Default to 100k reads for faster benchmarking

if [ -z "$GENOME" ] || [ -z "$READS" ]; then
    echo "Usage: $0 <genome.fa> <reads.fastq> [num_reads]"
    exit 1
fi

# Create benchmark output directory
BENCH_DIR="./benchmark-results-$(date +%Y%m%d-%H%M%S)"
mkdir -p "$BENCH_DIR"

# Create a subset of reads for faster testing
SUBSET_READS="$BENCH_DIR/reads_subset.fastq"
echo "Creating subset of $NUM_READS reads..."
head -n $((NUM_READS * 4)) "$READS" > "$SUBSET_READS"

# Results file
RESULTS="$BENCH_DIR/results.csv"
echo "level,index_size_bytes,index_time_sec,map_peak_rss_bytes,map_time_sec,mapped_reads" > "$RESULTS"

# Build fresh
echo "Building indexer and mapper..."
make clean && make -j4 2>/dev/null

# Benchmark each optimization level
for LEVEL in 0 1 2 3 4 5; do
    echo ""
    echo "========================================"
    echo "Testing optimization level $LEVEL"
    echo "========================================"
    
    INDEX_FILE="$BENCH_DIR/genome_L${LEVEL}.fmidx"
    OUTPUT_FILE="$BENCH_DIR/output_L${LEVEL}.sam"
    
    # Build index
    echo "Building index (level $LEVEL)..."
    INDEX_START=$(date +%s.%N)
    ./indexer -R "$GENOME" -I "$INDEX_FILE" -L "$LEVEL" 2>&1 | tee "$BENCH_DIR/indexer_L${LEVEL}.log"
    INDEX_END=$(date +%s.%N)
    INDEX_TIME=$(echo "$INDEX_END - $INDEX_START" | bc)
    INDEX_SIZE=$(stat -f%z "$INDEX_FILE" 2>/dev/null || stat --printf="%s" "$INDEX_FILE")
    
    echo "Index size: $INDEX_SIZE bytes ($((INDEX_SIZE / 1024 / 1024)) MB)"
    echo "Index time: ${INDEX_TIME}s"
    
    # Run mapper
    echo "Running mapper (level $LEVEL)..."
    MAP_START=$(date +%s.%N)
    
    # Run mapper with time command to get peak RSS
    if [ "$LEVEL" -eq 0 ]; then
        # Level 0 has genome in index, no need for -R
        /usr/bin/time -l ./mapper -I "$INDEX_FILE" -i "$SUBSET_READS" -o "$OUTPUT_FILE" -k 1 -t 1 2>&1 | tee "$BENCH_DIR/mapper_L${LEVEL}.log"
    else
        # Other levels need streaming FASTA
        /usr/bin/time -l ./mapper -I "$INDEX_FILE" -i "$SUBSET_READS" -o "$OUTPUT_FILE" -k 1 -t 1 -R "$GENOME" 2>&1 | tee "$BENCH_DIR/mapper_L${LEVEL}.log"
    fi
    
    MAP_END=$(date +%s.%N)
    MAP_TIME=$(echo "$MAP_END - $MAP_START" | bc)
    
    # Extract peak RSS from time output (macOS format)
    PEAK_RSS=$(grep "maximum resident set size" "$BENCH_DIR/mapper_L${LEVEL}.log" | awk '{print $1}')
    if [ -z "$PEAK_RSS" ]; then
        # Try to get from mapper's own output
        PEAK_RSS=$(grep "peak RSS:" "$BENCH_DIR/mapper_L${LEVEL}.log" | grep -oE '[0-9.]+ MiB' | head -1 | awk '{print int($1 * 1024 * 1024)}')
    fi
    
    # Count mapped reads
    MAPPED=$(grep -v "^\*" "$OUTPUT_FILE" 2>/dev/null | wc -l | tr -d ' ')
    
    echo "Peak RSS: $PEAK_RSS bytes ($((PEAK_RSS / 1024 / 1024)) MB)"
    echo "Map time: ${MAP_TIME}s"
    echo "Mapped reads: $MAPPED"
    
    # Record results
    echo "$LEVEL,$INDEX_SIZE,$INDEX_TIME,$PEAK_RSS,$MAP_TIME,$MAPPED" >> "$RESULTS"
done

echo ""
echo "========================================"
echo "BENCHMARK SUMMARY"
echo "========================================"
echo ""
cat "$RESULTS" | column -t -s ','
echo ""
echo "Results saved to: $RESULTS"
echo ""

# Generate comparison table
echo "Memory Savings Comparison (vs Level 0 baseline):"
echo "================================================"
BASELINE_RSS=$(grep "^0," "$RESULTS" | cut -d',' -f4)
BASELINE_IDX=$(grep "^0," "$RESULTS" | cut -d',' -f2)

for LEVEL in 1 2 3 4 5; do
    LINE=$(grep "^$LEVEL," "$RESULTS")
    if [ -n "$LINE" ]; then
        IDX_SIZE=$(echo "$LINE" | cut -d',' -f2)
        RSS=$(echo "$LINE" | cut -d',' -f4)
        
        IDX_SAVINGS=$(echo "scale=1; (1 - $IDX_SIZE / $BASELINE_IDX) * 100" | bc)
        RSS_SAVINGS=$(echo "scale=1; (1 - $RSS / $BASELINE_RSS) * 100" | bc)
        
        echo "Level $LEVEL: Index ${IDX_SAVINGS}% smaller, RSS ${RSS_SAVINGS}% lower"
    fi
done

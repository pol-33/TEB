# Mapper Progress Monitoring - Summary

## What Was Changed

I've added progress messages to the mapper to help you see what it's doing during execution.

### Changes Made:

1. **mapper.cpp** - Added progress messages at multiple levels:
   - Startup messages showing FM-index loading and configuration
   - Regular progress updates every 10,000 reads (was 100,000)
   - DEBUG-only verbose messages showing:
     - Every 1,000 reads being processed
     - Which chromosome is being searched (e.g., "chromosome 156/194")
     - Whether exact or inexact search is being used
     - How many search results were found

2. **Makefile** - Added DEBUG build support:
   - Use `make DEBUG=1` to build with verbose debug messages
   - Use `make` (without DEBUG=1) for production build with less verbose progress

3. **watch-progress.sh** - Created a monitoring script to watch progress in real-time

## How to Use

### Option 1: Non-DEBUG Build (Recommended for benchmarking)
```bash
make clean && make
INDEX=./genome.fmidx OUT_DIR=./bench-final BENCH_READS=0 RUN_BWA_CHECK=0 ./bench.sh
```

This will show progress every 10,000 reads without the detailed per-chromosome messages.

### Option 2: DEBUG Build (For detailed progress)
```bash
make clean && make DEBUG=1
INDEX=./genome.fmidx OUT_DIR=./bench-final BENCH_READS=0 RUN_BWA_CHECK=0 ./bench.sh
```

This shows detailed progress including:
- Every 1,000 reads
- Current chromosome being searched
- Search results for each chromosome

### Option 3: Watch Progress in Real-Time
In a separate terminal, run:
```bash
./watch-progress.sh
```

This will show a live updating view of the mapper's progress, refreshing every 2 seconds.

## What You'll See

### Non-DEBUG mode:
```
[mapper] loading FM-index from ./genome.fmidx
[mapper] index loaded, 194 chromosomes
[mapper] opening reads from ./bench-final/bench_reads.fastq
[mapper] starting read processing with k=1
[mapper] processed 10000 reads
[mapper] processed 20000 reads
...
[mapper] processed 1000000 reads total
[mapper] current RSS: XXX MiB, peak RSS: YYY MiB
```

### DEBUG mode (much more verbose):
```
[mapper] loading FM-index from ./genome.fmidx
[mapper] index loaded, 194 chromosomes
[mapper] opening reads from ./bench-final/bench_reads.fastq
[mapper] starting read processing with k=1
[mapper] processing read 0: ERR194147.787265110
[mapper-debug] read 0 forward search
[mapper-debug] searching chromosome 0/194
[mapper-debug] no exact match, starting inexact search
[mapper-debug] inexact search found 0 results
[mapper-debug] searching chromosome 1/194
...
```

## Note on Performance

The DEBUG build adds overhead due to the frequent I/O for progress messages. For accurate benchmarking, use the non-DEBUG build. The DEBUG build is primarily for troubleshooting and verifying the mapper is making progress.

## Current Status

Your benchmark is currently running with DEBUG enabled. You can monitor it by:
1. Checking the log file: `tail -f ./bench-final/mapper.stderr.log`
2. Using the watch script: `./watch-progress.sh`
3. Checking specific progress: `grep "processed.*reads" ./bench-final/mapper.stderr.log`

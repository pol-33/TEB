#!/bin/bash
# Watch mapper progress in real-time
# Usage: ./watch-progress.sh [logfile]

LOGFILE="${1:-./bench-final/mapper.stderr.log}"

if [[ ! -f "$LOGFILE" ]]; then
    echo "Log file not found: $LOGFILE"
    exit 1
fi

echo "Watching mapper progress from: $LOGFILE"
echo "Press Ctrl+C to stop"
echo ""

while true; do
    clear
    echo "=============== MAPPER PROGRESS ==============="
    echo ""
    
    # Show latest read progress
    echo "Latest read processed:"
    grep "processing read [0-9]" "$LOGFILE" | tail -1
    echo ""
    
    # Show latest 10k milestone
    echo "Latest 10k milestone:"
    grep "processed [0-9]* reads" "$LOGFILE" | tail -1
    echo ""
    
    # Show current chromosome being searched
    echo "Current chromosome search:"
    grep "searching chromosome" "$LOGFILE" | tail -1
    echo ""
    
    # Show last few messages
    echo "Last 5 log messages:"
    tail -5 "$LOGFILE"
    echo ""
    echo "==============================================="
    
    sleep 2
done

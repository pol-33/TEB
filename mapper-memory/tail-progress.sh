#!/bin/bash
# Tail the progress log to see updates in real-time
# Usage: ./tail-progress.sh [logfile]

LOGFILE="${1:-./benchmark-dataset-complet/mapper.stderr.log}"

if [[ ! -f "$LOGFILE" ]]; then
    echo "Waiting for log file to appear: $LOGFILE"
    while [[ ! -f "$LOGFILE" ]]; do
        sleep 1
    done
fi

echo "Showing real-time progress from: $LOGFILE"
echo "Press Ctrl+C to stop"
echo ""

# Show last 20 lines first to give context
tail -20 "$LOGFILE"
echo ""
echo "--- Following new output ---"
echo ""

# Now follow the file
tail -f "$LOGFILE" | grep --line-buffered -E "(processing read [0-9]|processed [0-9]+ reads|searching chromosome|mapper\])"

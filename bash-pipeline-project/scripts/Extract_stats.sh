#!/bin/bash

# Load configuration variables
source "$PROJECT_ROOT/config/config.env"
# Create results directory for demultiplexing stats (if needed)
mkdir -p "$DEMULTIPLEX_DIR"

# Ensure write permissions on the results directory
chmod u+w "$DEMULTIPLEX_DIR"

# Find the detailed log file (it has a pattern process_radtags.*.log)
LOGFILE=$(find "$DEMULTIPLEX_DIR" -name "process_radtags.*.log" -type f | head -1)

# Fallback to the generic log if detailed log not found
if [ -z "$LOGFILE" ] || [ ! -f "$LOGFILE" ]; then
    LOGFILE="${DEMULTIPLEX_DIR}/process_radtags.log"
fi

OUTFILE="${DEMULTIPLEX_DIR}/stats_per_sample.csv"

# Extract and format the data section from the log file
awk '
  BEGIN {FS="[ \t]+"}  # Use spaces or tabs as field separators
  /BEGIN per_barcode_raw_read_counts/ {capture=1; next}   # Start capturing data
  /END per_barcode_raw_read_counts/ {capture=0}           # Stop capturing data
  capture && $1 != "Barcode" {
    # Remove percent signs for numerical consistency
    gsub("%", "", $8)
    gsub("%", "", $10)
    gsub("%", "", $12)

    # Output fields as CSV
    print $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $7 "," $8 "," $9 "," $10 "," $11 "," $12
  }
' "$LOGFILE" >> "$OUTFILE"
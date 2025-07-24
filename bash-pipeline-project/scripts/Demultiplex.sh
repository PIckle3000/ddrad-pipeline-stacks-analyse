#!/bin/bash -e

#SBATCH --account=uoa02626
#SBATCH --job-name=Demultiplex
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --hint=nomultithread
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clem.cha87@gmail.com
#SBATCH -N 1


   module load Stacks/2.67-GCC-12.3.0

# Load configuration variables
source "$PROJECT_ROOT/config/config.env"

# Create results directory for demultiplexing
mkdir -p "$DEMULTIPLEX_DIR"

INPUT_DIR=$(dirname "$FASTQ1")
BARCODE_FILE="$BARCODE_PATH"
OUTPUT_DIR="$DEMULTIPLEX_DIR"
LOGFILE="${OUTPUT_DIR}/process_radtags.log"
CSV_STATS="${OUTPUT_DIR}/stats_per_sample.csv"

# Backup old log if it exists
if [ -f "$LOGFILE" ]; then
    mv "$LOGFILE" "${LOGFILE}.$(date +%Y%m%d_%H%M%S)"
fi

cd "$OUTPUT_DIR"
# Demultiplex raw reads
process_radtags -P -p "$INPUT_DIR" \
    -b "$BARCODE_FILE" \
    -o "$OUTPUT_DIR" \
    -c -q -r --inline-index --renz-1 hindIII --renz-2 mspI -D -t 100 -i gzfastq -y gzfastq \
   > "$LOGFILE" 2>&1

touch "$DEMULTIPLEX_DIR/process_radtags.done"
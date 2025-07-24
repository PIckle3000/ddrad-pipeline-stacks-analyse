#!/bin/bash -e

#SBATCH --account=uoa02626
#SBATCH --job-name=Alignment
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clem.cha87@gmail.com

# Load required modules
module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.21-GCC-12.3.0

# Load configuration variables
source "$PROJECT_ROOT/config/config.env"    

# Create results directory for alignment
mkdir -p "$ALIGN_DIR"

GENOME="$GENOME_PATH"
READS_DIR="$DEMULTIPLEX_DIR"

# Get sample names from demultiplexed files
samples=$(ls ${READS_DIR}/*.1.fq.gz | grep -v '\.rem\.1\.fq\.gz' | sed 's#.*/##' | sed 's/\.1\.fq\.gz//')


for s in $samples
do
    # Skip if alignment already exists
    if [ -f "$ALIGN_DIR/${s}_aln.sam" ]; then
        echo "Sample $s already aligned, skipping."
        continue
    fi

    echo "Processing sample: $s"

    # Align paired-end reads (double threads: -t 8)
    bwa mem -t 8 "$GENOME" \
        "$READS_DIR/${s}.1.fq.gz" "$READS_DIR/${s}.2.fq.gz" \
        > "$ALIGN_DIR/${s}_PEaln.sam" &

    # Align unpaired remnant reads 1 (double threads: -t 4)
    if [ -f "$READS_DIR/${s}.rem.1.fq.gz" ]; then
        bwa mem -t 4 "$GENOME" "$READS_DIR/${s}.rem.1.fq.gz" > "$ALIGN_DIR/${s}_rem1aln.sam" &
    fi
    # Align unpaired remnant reads 2 (double threads: -t 4)
    if [ -f "$READS_DIR/${s}.rem.2.fq.gz" ]; then
        bwa mem -t 4 "$GENOME" "$READS_DIR/${s}.rem.2.fq.gz" > "$ALIGN_DIR/${s}_rem2aln.sam" &
    fi

    wait #wait alignments to finish

    #Combine SAM files with only the first header and all alignments
    awk '/^@/ { print } !/^@/ { exit }' "$ALIGN_DIR/${s}_PEaln.sam" > "$ALIGN_DIR/${s}_aln.sam"
    for file in "$ALIGN_DIR/${s}_PEaln.sam" "$ALIGN_DIR/${s}_rem1aln.sam" "$ALIGN_DIR/${s}_rem2aln.sam"
    do
        [ -f "$file" ] && awk '!/^@/' "$file" >> "$ALIGN_DIR/${s}_aln.sam"
    done

done
touch "$ALIGN_DIR/alignment.done"
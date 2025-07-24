#!/bin/bash
# filepath: scripts/stacks.sh

#SBATCH --job-name=stacks_analysis
#SBATCH --output=stacks_analysis_%j.out
#SBATCH --error=stacks_analysis_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=5:00:00
#SBATCH --partition=genoa

source "$PROJECT_ROOT/config/config.env"
module load Stacks/2.67-GCC-12.3.0

#in and out directories
if [ -n "$USER_POPMAP" ] && [ -f "$USER_POPMAP" ]; then
    POPMAP="$USER_POPMAP"
else
    POPMAP="$PROJECT_ROOT/config/popmapdefault.txt"
fi
BAM_DIR="$PROJECT_ROOT/results/bam"
CATALOGUED_DIR="$PROJECT_ROOT/results/stacks/catalogued"
POPS_DIR="$PROJECT_ROOT/results/stacks/pops"

mkdir -p "$CATALOGUED_DIR"
mkdir -p "$POPS_DIR"

# 1. Analyse all samples with gstacks
gstacks -M "$POPMAP" \
        -I "$BAM_DIR" \
        -O "$CATALOGUED_DIR" \
        --min-mapq 20 --var-alpha 0.01 --gt-alpha 0.01 --max-clipped 0.1 -t 6

# 2. Run populations to get all loci for initial QC
populations -P "$CATALOGUED_DIR" \
            -M "$POPMAP" \
            -O "$POPS_DIR" \
            -t 6 --vcf --hwe --genepop --structure

touch "$STACKS_DIR/stacks.done"
#!/bin/bash
#SBATCH --job-name=populations_analysis
#SBATCH --output=populations_analysis_%j.out
#SBATCH --error=populations_analysis_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --partition=genoa

source "$PROJECT_ROOT/config/config.env"
module load Stacks/2.67-GCC-12.3.0

# Utilise la popmap passée par le pipeline (USER_POPMAP ou générée)
POPMAP="$POPMAP"

# Utilise le VCF filtré comme entrée
FILTERED_VCF="$PROJECT_ROOT/results/stacks/filtered/out.5.recode.vcf"
FILTERED_DIR="$PROJECT_ROOT/results/stacks/filtered"

mkdir -p "$FILTERED_DIR"

# Rerun populations sur le VCF filtré, single SNP, whitelist si besoin
populations -V "$FILTERED_VCF" \
            -M "$POPMAP" \
            -O "$FILTERED_DIR" \
            -p 2 -r 0.5 --write-single-snp \
            --vcf --hwe --genepop --structure \
            -t 2
# Copy/rename the main VCF output for downstream analysis
if [ -f "$FILTERED_DIR/out.5.recode.p.snps.vcf" ]; then
    cp "$FILTERED_DIR/out.5.recode.p.snps.vcf" "$FILTERED_DIR/populations.snps.vcf"
    touch "$FILTERED_DIR/populations.done"
else
    echo "ERROR: populations did not finish successfully or output file not found." >&2
    exit 1
fi
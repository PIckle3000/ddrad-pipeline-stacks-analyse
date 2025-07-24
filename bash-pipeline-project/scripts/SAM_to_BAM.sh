#!/bin/bash -e
#SBATCH --account=uoa02626
#SBATCH --job-name=Convert_SAM_to_BAM
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --hint=nomultithread
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clem.cha87@gmail.com

module load SAMtools/1.21-GCC-12.3.0

# Load config
source "$PROJECT_ROOT/config/config.env"

samdir="$ALIGN_DIR"
outbam="$BAM_DIR"

mkdir -p "$outbam"

for samfile in "$samdir"/*.sam; do
  sample=$(basename "$samfile" .sam)
  echo "Processing $sample..."

  # Convert SAM → BAM, then sort
  samtools view -@ 8 -h -S -b "$samfile" -o "${outbam}/${sample}.bam"

  samtools sort -@ 8 -m 4G "${outbam}/${sample}.bam" \
    -T "${outbam}/${sample}_tmp" \
    -o "${outbam}/${sample}_sorted.bam"

  # Index the sorted BAM
  samtools index "${outbam}/${sample}_sorted.bam"

  # Cleanup
  rm -f "${outbam}/${sample}_tmp".*

  echo "$sample done."
done
touch "$BAM_DIR/sam_to_bam.done"
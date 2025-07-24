#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --account=uoa02626
#SBATCH --partition=genoa
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=clem.cha87@gmail.com

#charge bwa module
module load BWA/0.7.18-GCC-12.3.0

# Load configuration variables
source "$PROJECT_ROOT/config/config.env"

# Go to the genome directory
cd "$(dirname "$GENOME_PATH")"
# Index the reference genome
bwa index "$(basename "$GENOME_PATH")"

touch "$GENOME_PATH.bwt.done"
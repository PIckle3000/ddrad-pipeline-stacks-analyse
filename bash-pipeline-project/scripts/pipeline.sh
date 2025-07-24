#!/bin/bash
#SBATCH --job-name=ddRAD_pipeline
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=26:00:00
#SBATCH --mem=512M
#SBATCH --partition=genoa
#SBATCH --account=uoa02626
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=clem.cha87@gmail.com

export PROJECT_ROOT="/nesi/nobackup/uoa02626/Zyphius_demu/bash-pipeline-project"
source "$PROJECT_ROOT/config/config.env"

# Create logs directory if it doesn't exist
mkdir -p "$PROJECT_ROOT/logs"

# Redirect all output to a global pipeline log (and keep SLURM logs)
PIPELINE_LOG="$PROJECT_ROOT/logs/pipeline_global.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# 1. BWA index
if [ ! -f "$GENOME_PATH.bwt.done" ]; then
    echo "$(timestamp) - BWA index not found, running indexing step..."
    sbatch --output="$PROJECT_ROOT/logs/index_genome_%j.out" \
           --error="$PROJECT_ROOT/logs/index_genome_%j.err" \
           "$PROJECT_ROOT/scripts/index_genome.sh"
    while [ ! -f "$GENOME_PATH.bwt.done" ]; do
    echo "$(timestamp) - Waiting for BWA index to finish..."
        sleep 60
    done
    echo "$(timestamp) - BWA index completed"
else
    echo "$(timestamp) - BWA index already exists, skipping indexing step."
fi

# 2. Demultiplexing
if [ ! -f "$DEMULTIPLEX_DIR/process_radtags.done" ]; then
    echo "$(timestamp) - Demultiplexing not done, running demultiplexing step..."
    sbatch --output="$PROJECT_ROOT/logs/Demultiplex_%j.out" \
           --error="$PROJECT_ROOT/logs/Demultiplex_%j.err" \
           "$PROJECT_ROOT/scripts/Demultiplex.sh"
    while [ ! -f "$DEMULTIPLEX_DIR/process_radtags.done" ]; do
        echo "$(timestamp) - Waiting for demultiplexing to finish..."
        sleep 60
    done
    echo "$(timestamp) - Demultiplexing completed"

    # 3. Extract stats (no wait)
    echo "$(timestamp) - Running stats extraction step..."
    sbatch "$PROJECT_ROOT/scripts/Extract_stats.sh"
else
    echo "$(timestamp) - Demultiplexing already done, skipping step."
fi



# 4. Alignment
if [ ! -f "$ALIGN_DIR/alignment.done" ]; then
    echo "$(timestamp) - Alignment not done, running alignment step..."
    sbatch --output="$PROJECT_ROOT/logs/align_samples_%j.out" \
           --error="$PROJECT_ROOT/logs/align_samples_%j.err" \
           "$PROJECT_ROOT/scripts/align_samples.sh"
    while [ ! -f "$ALIGN_DIR/alignment.done" ]; do
        echo "$(timestamp) - Waiting for alignment to finish..."
        sleep 60
    done
    echo "$(timestamp) - Alignment completed"
else
    echo "$(timestamp) - Alignment already done, skipping step."
fi

# 5. SAM to BAM
if [ ! -f "$BAM_DIR/sam_to_bam.done" ]; then
    echo "$(timestamp) - SAM to BAM conversion not done, running conversion step..."
    sbatch --output="$PROJECT_ROOT/logs/SAM_to_BAM_%j.out" \
           --error="$PROJECT_ROOT/logs/SAM_to_BAM_%j.err" \
           "$PROJECT_ROOT/scripts/SAM_to_BAM.sh"
    while [ ! -f "$BAM_DIR/sam_to_bam.done" ]; do
        echo "$(timestamp) - Waiting for SAM to BAM conversion to finish..."
        sleep 60
    done
    echo "$(timestamp) - SAM to BAM conversion completed"
else
    echo "$(timestamp) - SAM to BAM conversion already done, skipping step."
fi


# 6. Popmap
if [ -n "$USER_POPMAP" ] && [ -f "$USER_POPMAP" ]; then
    echo "$(timestamp) - User popmap provided: $USER_POPMAP"
    cp "$USER_POPMAP" "$POPMAP"
    touch "$POPMAP.done"
elif [ ! -f "$POPMAP.done" ]; then
    echo "$(timestamp) - Generating popmap..."
    chmod +x "$PROJECT_ROOT/scripts/generate_popmap_bai.sh"
    "$PROJECT_ROOT/scripts/generate_popmap_bai.sh"
    while [ ! -f "$POPMAP.done" ]; do
        echo "$(timestamp) - Waiting for popmap generation to finish..."
        sleep 10
    done
    echo "$(timestamp) - Popmap generated successfully"
else
    echo "$(timestamp) - Popmap already generated, skipping step."
fi

# 7. Stacks
if [ ! -f "$STACKS_DIR/stacks.done" ]; then
    echo "$(timestamp) - Stacks not run yet, running stacks step..."
    sbatch --output="$PROJECT_ROOT/logs/stacks_%j.out" \
           --error="$PROJECT_ROOT/logs/stacks_%j.err" \
           "$PROJECT_ROOT/scripts/stacks.sh"
    while [ ! -f "$STACKS_DIR/stacks.done" ]; do
        echo "$(timestamp) - Waiting for stacks analysis to finish..."
        sleep 120
    done
    echo "$(timestamp) - Stacks analysis completed"
else
    echo "$(timestamp) - Stacks already run, skipping step."
fi

echo "$(timestamp) - Pipeline completed."
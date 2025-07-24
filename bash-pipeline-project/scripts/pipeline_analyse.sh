#!/bin/bash
#SBATCH --job-name=ddrad_complete_workflow
#SBATCH --account=uoa02626
#SBATCH --time=6:00:00
#SBATCH --mem=25GB
#SBATCH --cpus-per-task=1
#SBATCH --partition=genoa
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err

# ============================================================================== 
# Set PROJECT_ROOT (user must set this before running, or edit here)
# ==============================================================================

# Option 1: User exports PROJECT_ROOT before sbatch
# Option 2: Set it here (edit this path for your project location)
export PROJECT_ROOT="/nesi/nobackup/uoa02626/Zyphius_demu/Analyse/BaseTest/bash-pipeline-project"

# ============================================================================== 
# Load configuration and modules
# ==============================================================================

CONFIG_PATH="$PROJECT_ROOT/config/config.env"
if [ -f "$CONFIG_PATH" ]; then
    source "$CONFIG_PATH"
    echo "Loaded configuration from $CONFIG_PATH"
else
    echo "ERROR: config.env not found at $CONFIG_PATH"
    exit 1
fi

module purge
module load R/4.3.2-foss-2023a
module load VCFtools

timestamp() { date '+%Y-%m-%d %H:%M:%S'; }

# ============================================================================== 
# Set up variables from config
# ==============================================================================

STACKS_DIR="${STACKS_DIR:-$PROJECT_ROOT/results/stacks/pops}"
FILTERED_DIR="$PROJECT_ROOT/results/stacks/filtered"
ANALYSIS_DIR="$PROJECT_ROOT/results/stacks/analysis"
SPECIES="$SPECIES"

# Popmap selection logic
if [ -n "$USER_POPMAP" ]; then
    POPMAP="$USER_POPMAP"
    echo "Using USER_POPMAP: $POPMAP"
else
    POPMAP="$POPMAP"
    echo "Using default POPMAP: $POPMAP"
fi

export USER_POPMAP
export POPMAP
export STACKS_DIR="$STACKS_DIR"
export ANALYSIS_DIR="$ANALYSIS_DIR"
export SPECIES="$SPECIES"
# ============================================================================== 
# Logging setup
# ==============================================================================

mkdir -p "$PROJECT_ROOT/logs"
LOG_DIR="$PROJECT_ROOT/logs"
SCRIPT_LOG="$LOG_DIR/pipeline_analyse_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$SCRIPT_LOG") 2>&1

echo "$(timestamp) - Starting complete ddRAD workflow"
echo "Project root: $PROJECT_ROOT"
echo "Stacks directory: $STACKS_DIR"
echo "Species: $SPECIES"
echo "Popmap: $POPMAP"
echo "Working directory: $(pwd)"


# ==============================================================================
# STEP 1: FILTERING
# ==============================================================================

echo "$(timestamp) - Step 1: VCF Filtering and QC..."
FILTERING_PREREQUISITES_OK=true

if [ ! -f "$STACKS_DIR/populations.snps.vcf" ]; then
    echo "$(timestamp) - ERROR: Stacks output not found at $STACKS_DIR/populations.snps.vcf"
    FILTERING_PREREQUISITES_OK=false
fi

if [ -f "$STACKS_DIR/populations.snps.vcf" ]; then
    VCF_SIZE=$(stat -c%s "$STACKS_DIR/populations.snps.vcf" 2>/dev/null || echo "0")
    if [ "$VCF_SIZE" -lt 1000 ]; then
        echo "$(timestamp) - ERROR: VCF file appears empty or corrupted (size: $VCF_SIZE bytes)"
        FILTERING_PREREQUISITES_OK=false
    else
        echo "$(timestamp) - VCF file size: $VCF_SIZE bytes (OK)"
    fi
fi

if [ "$FILTERING_PREREQUISITES_OK" = false ]; then
    echo "$(timestamp) - FATAL: Filtering prerequisites not met. Exiting."
    exit 1
fi

if [ -f "$FILTERED_DIR/filtering.done" ]; then
    echo "$(timestamp) - Filtering already completed."
fi

if [ ! -f "$FILTERED_DIR/filtering.done" ]; then
    echo "$(timestamp) - Starting VCF filtering and QC..."
    cd "$PROJECT_ROOT/scripts"
    Rscript filter_vcf_ddrad.R
    if [ $? -eq 0 ]; then
        FILTERING_SUCCESS=true
        REQUIRED_FILES=(
            "$FILTERED_DIR/out.5.recode.vcf"
            "$FILTERED_DIR/TruesGlobalQC.wl"
            "$FILTERED_DIR/missingness_per_sample.initial.pdf"
        )
        for file in "${REQUIRED_FILES[@]}"; do
            if [ ! -f "$file" ]; then
                echo "$(timestamp) - ERROR: Required output file not found: $file"
                FILTERING_SUCCESS=false
            fi
        done
        if [ "$FILTERING_SUCCESS" = true ]; then
            FINAL_SNPS=$(grep -v "^#" "$FILTERED_DIR/out.5.recode.vcf" | wc -l 2>/dev/null || echo "0")
            FINAL_LOCI=$(wc -l < "$FILTERED_DIR/TruesGlobalQC.wl" 2>/dev/null || echo "0")
            echo "$(timestamp) - Filtering completed successfully!"
            echo "$(timestamp) - Final SNPs: $FINAL_SNPS"
            echo "$(timestamp) - Final loci: $FINAL_LOCI"
            touch "$FILTERED_DIR/filtering.done"
            echo "$(timestamp) - Results available in: $FILTERED_DIR"
        else
            echo "$(timestamp) - ERROR: Filtering completed but required output files missing"
            exit 1
        fi
    else
        echo "$(timestamp) - ERROR: R filtering script failed"
        exit 1
    fi
fi

if [ -f "$FILTERED_DIR/filtering.done" ] && [ ! -f "$FILTERED_DIR/populations.done" ]; then
    echo "$(timestamp) - Running populations on filtered VCF..."
    sbatch --output="$PROJECT_ROOT/logs/populations_%j.out" \
           --error="$PROJECT_ROOT/logs/populations_%j.err" \
           "$PROJECT_ROOT/scripts/populations.sh"
    while [ ! -f "$FILTERED_DIR/populations.done" ]; do
        echo "$(timestamp) - Waiting for populations analysis to finish..."
        sleep 60
    done
    echo "$(timestamp) - Populations analysis completed"
else
    if [ -f "$FILTERED_DIR/populations.done" ]; then
        echo "$(timestamp) - Populations analysis already done, skipping step."
    fi
fi

# ==============================================================================
# STEP 2: POPULATION GENETIC ANALYSIS
# ==============================================================================

echo "$(timestamp) - Step 2: Population Genetic Analysis..."

ANALYSIS_PREREQUISITES_OK=true

if [ ! -f "$FILTERED_DIR/filtering.done" ]; then
    echo "$(timestamp) - ERROR: Filtering not completed (missing filtering.done flag)"
    ANALYSIS_PREREQUISITES_OK=false
fi

POPS_VCF="$FILTERED_DIR/populations.snps.vcf"
if [ -f "$POPS_VCF" ]; then
    echo "$(timestamp) - Found populations VCF: $POPS_VCF"
    FINAL_VCF="$POPS_VCF"
    VCF_FOUND=true
else
    echo "$(timestamp) - ERROR: populations.snps.vcf not found in $FILTERED_DIR"
    VCF_FOUND=false
    ANALYSIS_PREREQUISITES_OK=false
fi

if [ "$VCF_FOUND" = true ]; then
    VCF_LINES=$(grep -v "^#" "$FINAL_VCF" | wc -l 2>/dev/null || echo "0")
    if [ "$VCF_LINES" -lt 100 ]; then
        echo "$(timestamp) - WARNING: VCF file has only $VCF_LINES SNPs (might be too low)"
    else
        echo "$(timestamp) - VCF file contains $VCF_LINES SNPs (OK)"
    fi
fi

if [ "$ANALYSIS_PREREQUISITES_OK" = false ]; then
    echo "$(timestamp) - FATAL: Analysis prerequisites not met. Exiting."
    exit 1
fi

if [ -f "$ANALYSIS_DIR/analysis.done" ]; then
    echo "$(timestamp) - Analysis already completed, checking results..."
    REQUIRED_ANALYSIS_FILES=(
        "$ANALYSIS_DIR/${SPECIES}_final.vcf"
        "$ANALYSIS_DIR/${SPECIES}_population_assignments.csv"
        "$ANALYSIS_DIR/${SPECIES}_analysis_summary.txt"
    )
    ANALYSIS_COMPLETE=true
    for file in "${REQUIRED_ANALYSIS_FILES[@]}"; do
        if [ ! -f "$file" ]; then
            echo "$(timestamp) - Missing analysis file: $file"
            ANALYSIS_COMPLETE=false
        fi
    done
    if [ "$ANALYSIS_COMPLETE" = true ]; then
        echo "$(timestamp) - Analysis already completed successfully. Skipping."
    else
        echo "$(timestamp) - Analysis done flag exists but results incomplete. Re-running..."
        rm -f "$ANALYSIS_DIR/analysis.done"
    fi
fi

if [ ! -f "$ANALYSIS_DIR/analysis.done" ]; then
    echo "$(timestamp) - Found filtered data, proceeding with analysis..."

    # Clear any R package installation locks
    R_LIB_PATH="/home/$(whoami)/R/foss-2023a/4.3"
    if [ -d "$R_LIB_PATH" ]; then
        LOCK_DIRS=$(find "$R_LIB_PATH" -name "00LOCK-*" -type d 2>/dev/null)
        if [ ! -z "$LOCK_DIRS" ]; then
            echo "$(timestamp) - Found installation locks, removing them:"
            for lock_dir in $LOCK_DIRS; do
                echo "  Removing: $lock_dir"
                rm -rf "$lock_dir"
            done
            echo "$(timestamp) - Installation locks cleared."
        else
            echo "$(timestamp) - No installation locks found."
        fi
    else
        echo "$(timestamp) - Creating R library directory: $R_LIB_PATH"
        mkdir -p "$R_LIB_PATH"
    fi
    
    # Export FILTERED_VCF for R scripts
    export FILTERED_VCF="$FINAL_VCF"

    # Run modular R scripts for analysis
    echo "$(timestamp) - Running Genetic diversity..."
    Rscript GD.R --species "$SPECIES" --popmap "$POPMAP" --vcf "$FINAL_VCF" --outdir "$ANALYSIS_DIR"

    echo "$(timestamp) - Running FST analysis..."
    Rscript FST.r --species "$SPECIES" --popmap "$POPMAP" --vcf "$FINAL_VCF" --outdir "$ANALYSIS_DIR"

    echo "$(timestamp) - Running LEA analysis..."
    Rscript LEA.r --species "$SPECIES" --popmap "$POPMAP" --vcf "$FINAL_VCF" --outdir "$ANALYSIS_DIR"

    echo "$(timestamp) - Running private alleles analysis..."
    Rscript UL.r --species "$SPECIES" --popmap "$POPMAP" --vcf "$FINAL_VCF" --outdir "$ANALYSIS_DIR"

    # Check for essential output files
    REQUIRED_FILES=(
        "$ANALYSIS_DIR/${SPECIES}_final.vcf"
        "$ANALYSIS_DIR/${SPECIES}_population_assignments.csv"
        "$ANALYSIS_DIR/${SPECIES}_analysis_summary.txt"
    )
    ANALYSIS_SUCCESS=true
    for file in "${REQUIRED_FILES[@]}"; do
        if [ ! -f "$file" ]; then
            echo "$(timestamp) - ERROR: Required output file not found: $file"
            ANALYSIS_SUCCESS=false
        fi
    done
    if [ "$ANALYSIS_SUCCESS" = true ]; then
        FINAL_SAMPLES=$(tail -n +2 "$ANALYSIS_DIR/${SPECIES}_population_assignments.csv" | wc -l 2>/dev/null || echo "0")
        echo "$(timestamp) - ddRAD analysis completed successfully!"
        echo "$(timestamp) - Final samples analyzed: $FINAL_SAMPLES"
        touch "$ANALYSIS_DIR/analysis.done"
        echo "$(timestamp) - Results available in: $ANALYSIS_DIR"
    else
        echo "$(timestamp) - ERROR: Analysis completed but required output files missing"
        exit 1
    fi
fi

# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

echo ""
echo "$(timestamp) - ===== COMPLETE WORKFLOW FINISHED ====="
echo "$(timestamp) - All steps completed successfully!"
echo ""
echo "Results locations:"
echo "  - Filtered data: $FILTERED_DIR"
echo "  - Analysis results: $ANALYSIS_DIR"
echo ""
echo "Key output files:"
echo "  - ${SPECIES}_analysis_summary.txt: Overall summary"
echo "  - ${SPECIES}_population_assignments.csv: Population assignments"
echo "  - ${SPECIES}_structure_plot.pdf: Population structure visualization"

if [ -f "$ANALYSIS_DIR/${SPECIES}_analysis_summary.txt" ]; then
    echo ""
    echo "$(timestamp) - Analysis summary:"
    cat "$ANALYSIS_DIR/${SPECIES}_analysis_summary.txt"
fi

echo ""
echo "$(timestamp) - Pipeline completed successfully!"
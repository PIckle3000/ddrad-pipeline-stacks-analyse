# ddRAD Bash Pipeline

> **Developed by Clément CHAPTAL**

---

## Overview

This project provides a modular Bash pipeline for processing ddRAD sequencing data, from raw FASTQ files to population genetics-ready outputs.  
The pipeline is fully configurable via environment variables and automatically manages result directories for each step.

---

## Project Structure

```
bash-pipeline-project
├── scripts/
│   ├── pipeline.sh                # Main orchestrator script
│   ├── index_genome.sh            # BWA genome indexing
│   ├── Demultiplex.sh             # Demultiplexing with process_radtags
│   ├── Extract_stats.sh           # Extraction of demultiplexing stats
│   ├── align_samples.sh           # Read alignment with BWA
│   ├── SAM_to_BAM.sh              # SAM to BAM conversion and sorting
│   ├── fastq_pairs.sh             # Fastq-pair utility
│   ├── stacks.sh                  # Stacks analysis
│   ├── generate_popmap_bai.sh     # Popmap generation (auto or user-provided)
│   ├── complete_ddrad_analysis.R  # R script for population genetic analysis
│   ├── filter_vcf_ddrad.R         # R script for VCF filtering
│   └── pipeline_analyse.sh        # Bash script for launching R analysis
├── config/
│   ├── config.env                 # Main configuration file
│   └── popmapdefault.txt          # Default popmap (optional)
├── results/
│   ├── demultiplex/
│   │   ├── sample01.1.fq.gz            # R1 FASTQ for sample01 (demultiplexed, cleaned)
│   │   ├── sample01.2.fq.gz            # R2 FASTQ for sample01
│   │   ├── sample01.rem.1.fq.gz        # R1 reads removed (low quality/barcode) for sample01
│   │   ├── sample01.rem.2.fq.gz        # R2 reads removed for sample01
│   │   ├── process_radtags.done        # Flag file: demultiplexing completed
│   │   ├── process_radtags.log         # Log file from process_radtags
│   │   ├── process_radtags.<run>.log   # Log for this run/species
│   │   ├── stats_per_sample.csv        # Summary stats per sample (reads kept, removed, etc.)
│   │   └── ...
│   ├── align/
│   │   ├── sample01_aln.sam            # Alignment SAM file for sample01 merged
│   │   ├── sample01_PEaln.sam          # Alignment SAM file for sample01 paired-end
│   │   ├── sample01_rem1aln.sam        # Alignment SAM file for sample01, removed reads set 1
│   │   ├── sample01_rem2aln.sam        # Alignment SAM file for sample01, removed reads set 2
│   │   ├── alignment.done              # Flag file: alignment completed
│   │   └── ...
│   ├── bam/
│   │   ├── sample01_aln.bam                # BAM file for sample01 merged
│   │   ├── sample01_aln_sorted.bam         # Sorted BAM file for sample01
│   │   ├── sample01_aln_sorted.bam.bai     # BAM index for sample01
│   │   ├── sample01_PEaln.bam              # BAM file for sample01 paired-end
│   │   ├── sample01_PEaln_sorted.bam       # Sorted BAM file for sample01 paired-end
│   │   ├── sample01_PEaln_sorted.bam.bai   # BAM index for sample01 paired-end
│   │   ├── sample01_rem1aln.bam            # BAM file for sample01, removed reads set 1
│   │   ├── sample01_rem1aln_sorted.bam     # Sorted BAM file for sample01, removed reads set 1
│   │   ├── sample01_rem1aln_sorted.bam.bai # BAM index for sample01, removed reads set 1
│   │   ├── sample01_rem2aln.bam            # BAM file for sample01, removed reads set 2
│   │   ├── sample01_rem2aln_sorted.bam     # Sorted BAM file for sample01, removed reads set 2
│   │   ├── sample01_rem2aln_sorted.bam.bai # BAM index for sample01, removed reads set 2
│   │   ├── sam_to_bam.done                 # Flag file: BAM conversion completed
│   │   └── ...
│   └── stacks/             
│       ├── catalogued/
│       │   ├── catalog.calls           # Catalog calls file (locus/allele calls)
│       │   ├── catalog.chrs.tsv        # Chromosome assignments for catalog loci
│       │   ├── catalog.fa.gz           # FASTA of catalog loci
│       │   ├── gstacks.log             # Log file from gstacks
│       │   ├── gstacks.log.distribs    # Distribution log from gstacks
│       │   └── ...
│       ├── pops/
│       │   ├── populations.haplotypes.tsv        # Haplotypes table
│       │   ├── populations.haps.genepop          # Haplotypes in Genepop format
│       │   ├── populations.haps.vcf              # Haplotypes in VCF format
│       │   ├── populations.hapstats.tsv          # Haplotypes statistics
│       │   ├── populations.log                   # Log file from populations
│       │   ├── populations.log.distribs          # Distribution log from populations
│       │   ├── populations.snps.genepop          # SNPs in Genepop format
│       │   ├── populations.snps.vcf              # SNPs in VCF format
│       │   ├── populations.structure             # Structure format output
│       │   ├── populations.sumstats.tsv          # Summary statistics
│       │   ├── populations.sumstats_summary.tsv  # Summary statistics (summary)
│       │   ├── stacks.done  
│       │   └── ...                               # Flag file: populations step
│       ├── filtered/                             # Filtered VCF and QC/intermediate files
│       │   ├── populations.snps.vcf              # Final filtered VCF file for downstream analysis
│       │   ├── populations.structure             # STRUCTURE format file for population structure inference
│       │   ├── populations.snps.genepop          # SNPs in Genepop format for population genetics software
│       │   ├── populations.sumstats.tsv          # Per-locus summary statistics after filtering (tab-separated)
│       │   ├── populations.sumstats_summary.tsv  # Summary of per-locus statistics after filtering
│       │   ├── populations.haplotypes.tsv        # Table of filtered haplotypes for each sample
│       │   ├── populations.hapstats.tsv          # Haplotype statistics after filtering
│       │   ├── populations.log                   # Log file from the populations module (filtering and processing details)
│       │   ├── populations.log.distribs          # Distribution log from populations (e.g., depth, missingness)
│       │   ├── meandepth_per_locus.S5initial.pdf # PDF plot of mean sequencing depth per locus (step 5, initial)
│       │   ├── missingness_per_sample.initial.pdf# PDF plot of missing data per sample before filtering
│       │   ├── missingness_per_sample.step7.pdf  # PDF plot of missing data per sample after filtering
│       │   ├── out.*.imiss / out.*.ldepth        # Intermediate QC files (missing data, depth per sample/locus)
│       │   ├── remove.*.sites / remove.*.inds    # Lists of loci or individuals removed at each filtering step
│       │   ├── whitelist*.loci / *.wl            # Whitelist of loci retained after filtering
│       │   ├── bald.out.*.recode.vcf             # Filtered VCFs for specific subsets or test cases
│       │   ├── filtering.done                    # Flag file: filtering pipeline completed successfully
│       │   ├── populations.done                  # Flag file: populations filtering step completed
│       │   └── ...
│       ├── analysis/                             # Downstream population genetic analysis outputs
│       │   ├── <species>_final.vcf               # Filtered VCF used for analysis
│       │   ├── <species>_population_assignments.csv      # Population assignments after filtering
│       │   ├── <species>_heterozygosity_by_population.csv# Heterozygosity stats per population
│       │   ├── <species>_heterozygosity_tests.csv        # Heterozygosity test results (if 2 pops)
│       │   ├── <species>_analysis_summary.txt            # Analysis summary report
│       │   ├── fst_results.csv                           # Pairwise FST values
│       │   ├── <species>_private_alleles_filtered_subset.csv # Private alleles per population
│       │   ├── <species>_snmf_cross_entropy.pdf          # LEA sNMF cross-entropy plot
│       │   ├── <species>_clustering_results.csv          # LEA Q-matrix (if K > 1)
│       │   ├── <species>_structure_plot.pdf              # LEA structure plot (if K > 1)
│       │   ├── <species>_data.geno                       # LEA genotype file
│       │   ├── <species>_data.snmfProject / .snmf        # LEA project files
│       │   └── ...                                       # Other analysis outputs
│
├── logs/                          # Log files for pipeline and jobs
└── README.md                      # Project documentation
```

---

## Getting Started

### Prerequisites

- Bash shell (Linux/Unix)
- SLURM workload manager (for job submission)
- Required bioinformatics tools: BWA, Stacks, SAMtools

---

### Installation

1. **Clone the repository:**
   ```
   scp nesi
   cd bash-pipeline-project
   ```

2. **Edit the configuration file:**
   ```
   nano config/config.env
   ```
   - Set all paths to your genome, barcodes, FASTQ files, and output directories as needed.
   - If you want to provide your own popmap, set the `USER_POPMAP` variable to the path of your file. Leave it empty to auto-generate with one population (some analyse is not possible with only one population)

3. **Edit the pipeline root path:**
   - In `scripts/pipeline.sh` and `scripts/pipeline_analyse.sh`, set the `PROJECT_ROOT` variable at the top of the script to match the root directory of your project on your system.

---

### Running the Pipeline

Submit the main pipeline script via SLURM:
```
sbatch scripts/pipeline.sh
```
Or, for local testing (not recommended for full runs):
```
bash scripts/pipeline.sh
```

---

### Output

All results and intermediate files are stored in the `results/` directory, organized by step. Logs are saved in the `logs/` directory.

---

## Popmap Handling

- By default, the pipeline auto-generates a popmap file.
- To use your own popmap, set the `USER_POPMAP` variable in `config/config.env` to the path of your file.

---

## Modifying Variables

All pipeline behavior is controlled via `config/config.env` and the root is control by pipeline.sh
Edit the congig file to change input files, output locations, or pipeline options without modifying the scripts.

---

## Troubleshooting

- **Line endings:**  
  Ensure all scripts use Unix (LF) line endings. Convert with `dos2unix *.sh` if needed.
- **Logs:**  
  Check the `logs/` directory for detailed error and progress messages.
- **Modules:**  
  Make sure all required modules/tools are available on your system or cluster.

---

### Genetic Analysis Pipeline

After running the main bash pipeline, you can perform population genetic analyses using the R script:

To launch the genetic analysis, submit the following command from the root directory:
```
sbatch pipeline_analyse.sh
```
This will run the R analysis script as a SLURM job. The script will analyze the filtered VCF file and generate summary statistics, population assignments, FST, structure plots, and other results in the `results/stacks/analysis/` directory.

#### 4. Output
- Results are written in the `results/stacks/analysis/` directory (or as configured).
- Key output files: population assignments, diversity stats, FST, structure plots, summary report, etc.
- Some of the Analyse are on the `populations.sumstats_summary.tsv` check it for have complementary analysis.

#### 5. Notes
- If you have more than two populations, the script is robust to missing data and will include all populations in the summary tables.

---

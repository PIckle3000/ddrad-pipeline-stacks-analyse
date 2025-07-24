#!/usr/bin/env Rscript

# ddRAD Analysis Script (Analysis Only)
# Single-species population genetic analyses using pre-filtered data
# No filtering - uses output from separate filtering pipeline

# Load required libraries
# Note: dartR is not available on NeSI R 4.3.2
# Functions dependent on this package are commented out or have alternatives
cat("=== INITIALIZING PACKAGES ===\n")

# Define personal NeSI R library path
user_lib <- file.path("/home/acch270", "R", "foss-2023a", "4.3")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
  cat("Created personal R library at:", user_lib, "\n")
}
.libPaths(c(user_lib, .libPaths()))
cat(
  "Using library paths:\n",
  paste0(" - ", .libPaths(), collapse = "\n"),
  "\n\n"
)

# List of required packages
required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "reshape2",
  "vcfR", "adegenet", "poppr", "ape", "pegas",
  "hierfstat", "StAMPP", "LEA"
)

# Function to install missing packages in personal library
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(" - Package missing:", pkg, " → attempting to install...\n")
    tryCatch(
      {
        install.packages(pkg, repos = "https://cran.r-project.org", lib = user_lib)
        library(pkg, character.only = TRUE, lib.loc = user_lib)
        cat("   ✓", pkg, "installed and loaded from", user_lib, "\n")
      },
      error = function(e) {
        cat("   ✗ Installation error for", pkg, ":", e$message, "\n")
      }
    )
  } else {
    cat("   ✓", pkg, "already available and loaded\n")
  }
}

cat("Loading required packages:\n")
invisible(lapply(required_packages, install_if_missing))

cat("\n=== Packages initialized ===\n\n")
cat("\n=== ddRAD ANALYSIS PIPELINE ===\n")
cat("Starting population genetic analysis on pre-filtered data\n")

# Get species name from environment or use default
species_name <- Sys.getenv("SPECIES", "ziphius")
cat("Species:", species_name, "\n")

# Set up directories
if (!is.null(Sys.getenv("STACKS_DIR", unset = NA)) && !is.na(Sys.getenv("STACKS_DIR", unset = NA))) {
  stacks_base_dir <- dirname(Sys.getenv("STACKS_DIR"))
  filter_results_dir <- file.path(stacks_base_dir, "filtered")
  analysis_results_dir <- file.path(stacks_base_dir, "analysis")
} else {
  filter_results_dir <- "../results/stacks/filtered"
  analysis_results_dir <- "../results/stacks/analysis"
}

# Create analysis directory
if (!dir.exists(analysis_results_dir)) {
  dir.create(analysis_results_dir, recursive = TRUE)
  cat("Created analysis directory:", analysis_results_dir, "\n")
}

# Check if filtering results exist
if (!dir.exists(filter_results_dir)) {
  stop("ERROR: Filtering results directory not found. Run filtering pipeline first.")
}

# Use the final populations VCF file (after populations step)
final_vcf <- file.path(filter_results_dir, "populations.snps.vcf")
if (!file.exists(final_vcf)) {
  stop("ERROR: populations.snps.vcf not found in ", filter_results_dir)
}
# Copy the filtered VCF to the analysis directory with the expected name
analysis_vcf <- file.path(analysis_results_dir, paste0(species_name, "_final.vcf"))
file.copy(final_vcf, analysis_vcf, overwrite = TRUE)
cat("Copied filtered VCF to analysis directory:", analysis_vcf, "\n")

cat("Loading pre-filtered VCF data for population genetic analysis...\n")

# Load VCF with vcfR
vcf_data <- read.vcfR(analysis_vcf, verbose = FALSE)
cat("Loaded VCF with", nrow(vcf_data@fix), "SNPs and", ncol(vcf_data@gt) - 1, "samples\n")

# Convert to different formats for analysis
cat("Converting to different data formats...\n")
gl_data <- vcfR2genlight(vcf_data)

# Convert to genind format (using vcfR method)
cat("Converting to genind using vcfR\n")
gi_data <- vcfR2genind(vcf_data)

# Basic population assignment (user should modify based on their data)
n_samples_final <- ncol(vcf_data@gt) - 1
sample_names_final <- colnames(vcf_data@gt)[-1]

# --- Population assignment using popmap file from config ---
# Try to get popmap path from environment/config
get_popmap_path <- function() {
  # Try USER_POPMAP first
  user_popmap <- Sys.getenv("USER_POPMAP", unset = NA)
  if (!is.na(user_popmap) && user_popmap != "") {
    return(user_popmap)
  }
  # Try POPMAP from config.env
  popmap_env <- Sys.getenv("POPMAP", unset = NA)
  if (!is.na(popmap_env) && popmap_env != "") {
    return(popmap_env)
  }
  # Fallback: look for popmap in config directory
  config_dir <- file.path(dirname(dirname(getwd())), "config")
  popmap_default <- file.path(config_dir, "popmapdefault.txt")
  if (file.exists(popmap_default)) {
    return(popmap_default)
  }
  return(NA)
}

popmap_path <- get_popmap_path()
if (is.na(popmap_path) || !file.exists(popmap_path)) {
  stop("No valid popmap file found. Please set USER_POPMAP or POPMAP in config.env and rerun. Analysis aborted.")
}
cat("Using popmap file:", popmap_path, "\n")

# Read popmap file (tab or space separated)
popmap_df <- tryCatch(
  {
    read.table(popmap_path, header = FALSE, stringsAsFactors = FALSE)
  },
  error = function(e) {
    stop("Could not read popmap file: ", popmap_path)
  }
)
colnames(popmap_df) <- c("Sample", "Population")

# Match VCF samples to popmap
new_strata <- setNames(popmap_df$Population[match(sample_names_final, popmap_df$Sample)], sample_names_final)
if (any(is.na(new_strata))) {
  missing_samples <- sample_names_final[is.na(new_strata)]
  stop(paste0("The following samples are missing from the popmap file: ", paste(missing_samples, collapse = ", "), ". Analysis aborted."))
}

# Population assignment (direct approach)
cat("Starting population assignment (popmap)\n")
pop(gi_data) <- new_strata
hs_data <- genind2hierfstat(gi_data)
cat("Initial population assignment complete\n")

# Convert back to genlight format using adegenet
cat("Converting genind to genlight using adegenet\n")
if ("genind2genlight" %in% ls("package:adegenet")) {
  gl_data <- adegenet::genind2genlight(gi_data)
} else {
  cat("Warning: genind2genlight() not found in adegenet. Skipping conversion to genlight.\n")
  gl_data <- NULL
}

# Save population assignments
pop_assignments <- data.frame(
  Sample = sample_names_final,
  Population = new_strata
)

# Filter populations with at least 2 samples and exclude 'Unknown' immediately after population assignment
cat("Filtering valid populations (>=2 samples, excluding 'Unknown')\n")
tab <- table(new_strata)
valid_pops <- names(tab)[tab >= 2 & names(tab) != "Unknown"]
filtered_indices <- which(new_strata %in% valid_pops)
filtered_gi_data <- gi_data[filtered_indices]
filtered_strata <- new_strata[filtered_indices]
filtered_sample_names <- sample_names_final[filtered_indices]
cat("Population filtering complete\n")

# Use filtered_gi_data and filtered_strata everywhere below
# Update population assignments export
cat("Exporting filtered population assignments\n")
pop_assignments <- data.frame(
  Sample = filtered_sample_names,
  Population = filtered_strata
)
pop_assignments <- pop_assignments[!duplicated(pop_assignments$Sample), ]

# Use filtered_gi_data and filtered_strata for all downstream analyses
pop(filtered_gi_data) <- filtered_strata
hs_data <- genind2hierfstat(filtered_gi_data)
sample_names_final <- filtered_sample_names
new_strata <- filtered_strata

# ==============================================================================
# GENETIC DIVERSITY ANALYSES
# ==============================================================================

cat("\nCalculating genetic diversity statistics\n")
# Basic statistics
basic_stats <- basic.stats(hs_data, diploid = TRUE)
cat("Basic statistics calculated\n")

# Observed heterozygosity (Ho)
ho_stats <- basic_stats$Ho
ho_melt <- melt(ho_stats)[, 2:3]
colnames(ho_melt) <- c("pop", "ho")

# Expected heterozygosity (Hs)
hs_stats <- basic_stats$Hs
hs_melt <- melt(hs_stats)[, 2:3]
colnames(hs_melt) <- c("pop", "hs")

# Statistical tests for heterozygosity differences
if (length(unique(new_strata)) > 1) {
  # Get population names
  pop_names <- unique(new_strata)

  # Calculate means for each population
  ho_means <- aggregate(ho ~ pop, data = ho_melt, FUN = mean)
  hs_means <- aggregate(hs ~ pop, data = hs_melt, FUN = mean)

  cat("Observed Heterozygosity (Ho) by population:\n")
  print(ho_means)
  cat("Expected Heterozygosity (Hs) by population:\n")
  print(hs_means)

  # Statistical tests for Ho (following Mmirus methodology)
  if (length(pop_names) == 2) {
    ho_pop1 <- ho_melt[ho_melt$pop == pop_names[1], "ho"]
    ho_pop2 <- ho_melt[ho_melt$pop == pop_names[2], "ho"]

    # Kolmogorov-Smirnov test for Ho
    ho_ks_test <- ks.test(ho_pop1, ho_pop2)
    cat("\nKolmogorov-Smirnov test for Ho differences:\n")
    cat("D =", ho_ks_test$statistic, ", p-value =", ho_ks_test$p.value, "\n")

    # t-test for Ho
    ho_t_test <- t.test(ho_pop1, ho_pop2)
    cat("t-test for Ho differences:\n")
    cat("t =", ho_t_test$statistic, ", df =", ho_t_test$parameter, ", p-value =", ho_t_test$p.value, "\n")

    # Statistical tests for Hs
    hs_pop1 <- hs_melt[hs_melt$pop == pop_names[1], "hs"]
    hs_pop2 <- hs_melt[hs_melt$pop == pop_names[2], "hs"]

    # Kolmogorov-Smirnov test for Hs
    hs_ks_test <- ks.test(hs_pop1, hs_pop2)
    cat("\nKolmogorov-Smirnov test for Hs differences:\n")
    cat("D =", hs_ks_test$statistic, ", p-value =", hs_ks_test$p.value, "\n")

    # t-test for Hs
    hs_t_test <- t.test(hs_pop1, hs_pop2)
    cat("t-test for Hs differences:\n")
    cat("t =", hs_t_test$statistic, ", df =", hs_t_test$parameter, ", p-value =", hs_t_test$p.value, "\n")

    # Save statistical test results
    het_stats <- data.frame(
      Population = pop_names,
      Ho_mean = ho_means$ho,
      Hs_mean = hs_means$hs
    )

    het_test_results <- data.frame(
      Test = c("Ho_KS_test", "Ho_t_test", "Hs_KS_test", "Hs_t_test"),
      Statistic = c(ho_ks_test$statistic, ho_t_test$statistic, hs_ks_test$statistic, hs_t_test$statistic),
      P_value = c(ho_ks_test$p.value, ho_t_test$p.value, hs_ks_test$p.value, hs_t_test$p.value),
      DF = c(NA, ho_t_test$parameter, NA, hs_t_test$parameter)
    )

    #
  } else {
    cat("Statistical tests only implemented for 2 populations. Using basic summary statistics.\n")

    # For more than 2 populations, robustly merge Ho and Hs means by population
    # This ensures all populations are included, even if missing from one stat
    het_stats <- dplyr::full_join(ho_means, hs_means, by = "pop")
    colnames(het_stats) <- c("Population", "Ho_mean", "Hs_mean")
    write.csv(het_stats, paste0(species_name, "_heterozygosity_by_population.csv"), row.names = FALSE)
  }
}



# ==============================================================================
# SUMMARY REPORT
# ==============================================================================

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Species/Dataset:", species_name, "\n")
cat("Final number of samples:", n_samples_final, "\n")

cat("Populations:", paste(unique(new_strata), collapse = ", "), "\n")

cat("\nProduced files by this script:\n")
main_outputs <- c(
  paste0(species_name, "_analysis_summary.txt"),
  paste0(species_name, "_population_assignments.csv"),
  paste0(species_name, "_heterozygosity_by_population.csv"),
  paste0(species_name, "_heterozygosity_tests.csv")
)
for (f in main_outputs) {
  if (file.exists(f)) cat("  *", f, "\n")
}

cat("\nAdditional analyses (HWE, FST, clustering) will be produced by specialized scripts.\n")
cat("Check the analysis directory after pipeline completion for all results.\n")


# Export population assignments (filtered)
write.csv(pop_assignments, file.path(analysis_results_dir, paste0(species_name, "_population_assignments.csv")), row.names = FALSE)
cat("Saved population assignments:", file.path(analysis_results_dir, paste0(species_name, "_population_assignments.csv")), "\n")
cat("Export complete\n")

# Export heterozygosity statistics and test results
write.csv(het_stats, file.path(analysis_results_dir, paste0(species_name, "_heterozygosity_by_population.csv")), row.names = FALSE)
if (exists("het_test_results")) {
  write.csv(het_test_results, file.path(analysis_results_dir, paste0(species_name, "_heterozygosity_tests.csv")), row.names = FALSE)
}
cat("Saved heterozygosity statistics and test results\n")

# Export summary report
summary_file <- file.path(analysis_results_dir, paste0(species_name, "_analysis_summary.txt"))
sink(summary_file)
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Species/Dataset:", species_name, "\n")
cat("Final number of samples:", n_samples_final, "\n")
cat("Populations:", paste(unique(new_strata), collapse = ", "), "\n")
sink()
cat("Saved summary report to:", summary_file, "\n")

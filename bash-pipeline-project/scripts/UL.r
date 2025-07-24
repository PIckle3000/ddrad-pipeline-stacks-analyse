#!/usr/bin/env Rscript

library(vcfR)
library(adegenet)
library(poppr)

analysis_results_dir <- Sys.getenv("ANALYSIS_DIR")
species_name <- Sys.getenv("SPECIES")
popmap_path <- Sys.getenv("POPMAP")

cat("----- Starting imbalance loci analysis -----\n")
vcf_file <- file.path(analysis_results_dir, paste0(species_name, "_final.vcf"))
cat("Reading VCF file:", vcf_file, "\n")
vcf <- read.vcfR(vcf_file)
cat("VCF file loaded.\n")

cat("Converting VCF to genind object...\n")
gi <- vcfR2genind(vcf)
cat("Conversion done.\n")

cat("Reading population map:", popmap_path, "\n")
popmap_df <- read.table(popmap_path, header = FALSE, stringsAsFactors = FALSE)
colnames(popmap_df) <- c("Sample", "Population")
sample_names <- colnames(vcf@gt)[-1]
populations <- setNames(popmap_df$Population[match(sample_names, popmap_df$Sample)], sample_names)
if (any(is.na(populations))) stop("Samples missing in popmap!")
pop(gi) <- populations
cat("Population assignment completed.\n")

cat("Filtering populations with at least 2 individuals...\n")
pop_sizes <- table(pop(gi))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
keep_idx <- which(pop(gi) %in% valid_pops)
gi <- gi[keep_idx]
cat("Filtering done.\n")

cat("Individuals remaining after filtering:", nInd(gi), "\n")
cat("Loci remaining after filtering:", nLoc(gi), "\n")

if (nLoc(gi) == 0) {
    cat("No loci remaining after filtering. Analysis stopped.\n")
    quit(save = "no", status = 0)
}
if (!validObject(gi)) stop("genind object is invalid after filtering.")

cat("Population sizes after filtering:\n")
print(table(pop(gi)))

cat("Calculating private alleles...\n")
if (length(unique(pop(gi))) > 1 && nLoc(gi) > 0) {
    private_alleles_res <- private_alleles(gi, alleles ~ pop)
    write.csv(private_alleles_res, file.path(analysis_results_dir, paste0(species_name, "_private_alleles_filtered_subset.csv")), row.names = TRUE)
    cat("Private alleles results saved.\n")
} else {
    cat("Private alleles analysis skipped: not enough populations or loci.\n")
}

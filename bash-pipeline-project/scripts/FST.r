#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(vcfR)
    library(adegenet)
    library(hierfstat)
    library(dplyr)
})

cat("----- Starting FST analysis  -----\n")

# Use environment variables for modularity
vcf_path <- Sys.getenv("FILTERED_VCF", "")
popmap_path <- Sys.getenv("POPMAP", "")
analysis_results_dir <- Sys.getenv("ANALYSIS_DIR", ".")
output_path <- file.path(analysis_results_dir, "fst_results.csv")

cat("Reading VCF file:", vcf_path, "\n")
vcf <- read.vcfR(vcf_path, verbose = FALSE)
cat("VCF loaded with", nrow(vcf@fix), "variants and", ncol(vcf@gt) - 1, "individuals\n")

cat("Converting VCF to genind...\n")
gi <- vcfR2genind(vcf)
cat("Genind created with", nInd(gi), "individuals and", nLoc(gi), "loci\n")

cat("Reading popmap:", popmap_path, "\n")
popmap <- read.table(popmap_path, header = FALSE, sep = "\t", col.names = c("ind", "pop"))
popmap <- popmap[popmap$ind %in% indNames(gi), ]
gi <- gi[popmap$ind] # Reorder
pop(gi) <- as.factor(popmap$pop)
cat("Assigned populations:\n")
print(table(pop(gi)))

pop_sizes <- table(pop(gi))
valid_pops <- names(pop_sizes[pop_sizes >= 2])
gi <- gi[pop(gi) %in% valid_pops, ]
cat("Retained populations (>=2 individuals):\n")
print(table(pop(gi)))
cat("Number of individuals after filtering:", nInd(gi), "\n")
cat("Number of loci (SNPs):", nLoc(gi), "\n")

cat("Does gi contain NA?", any(is.na(tab(gi))), "\n")
gi_tab <- tab(gi)
gi_tab[is.na(gi_tab)] <- 0
gi <- new("genind",
    tab = gi_tab,
    pop = pop(gi),
    ind.names = indNames(gi),
    loc.names = locNames(gi),
    ploidy = ploidy(gi),
    type = gi@type
)
cat("Does gi contain NA after replacement?", any(is.na(tab(gi))), "\n")

cat("Converting to hierfstat format...\n")
gi_hierf <- genind2hierfstat(gi)
cat("Computing pairwise FST with hierfstat...\n")
fst_matrix <- pairwise.WCfst(gi_hierf)

cat("Saving results to:", output_path, "\n")
write.csv(as.data.frame(fst_matrix), output_path, row.names = TRUE)
cat("FST analysis complete.\n")

#!/usr/bin/env Rscript

library(vcfR)
library(LEA)
library(reshape2)
library(ggplot2)

analysis_results_dir <- Sys.getenv("ANALYSIS_DIR")
species_name <- Sys.getenv("SPECIES")
vcf_file <- file.path(analysis_results_dir, paste0(species_name, "_final.vcf"))
geno_file <- file.path(analysis_results_dir, paste0(species_name, "_data.geno"))
cat("----- Starting LEA analysis -----\n")
cat
cat(" Converting VCF to .geno file...\n")
vcf2geno(input.file = vcf_file, output.file = geno_file, force = TRUE)
cat(".geno file created:", geno_file, "\n")

cat(" Running sNMF clustering...\n")
obj <- snmf(geno_file,
    K = 1:4,
    repetitions = 10,
    project = "new",
    tolerance = 1e-5,
    entropy = TRUE,
    ploidy = 2
)
cat("sNMF clustering finished.\n")

cat(" Saving cross-entropy plot...\n")
pdf(file.path(analysis_results_dir, paste0(species_name, "_snmf_cross_entropy.pdf")))
plot(obj, cex = 1.2, col = "lightblue", pch = 19)
dev.off()
cat("Cross-entropy plot saved.\n")

cat(" Determining best K value...\n")
Kvec <- 1:4
ce <- sapply(Kvec, function(k) cross.entropy(obj, K = k))
best_k <- Kvec[which.min(ce)]
cat("Best K:", best_k, "\n")

if (!is.na(best_k) && best_k > 1) {
    cat(" Extracting Q-matrix and saving clustering results...\n")
    best_run <- which.min(cross.entropy(obj, K = best_k))
    qmatrix <- Q(obj, K = best_k, run = best_run)

    sample_names <- colnames(read.vcfR(vcf_file)@gt)[-1]

    clustering_res <- data.frame(Sample = sample_names, qmatrix)
    write.csv(clustering_res, file.path(analysis_results_dir, paste0(species_name, "_clustering_results.csv")), row.names = FALSE)
    cat("Clustering results saved.\n")

    cat(" Creating and saving structure plot...\n")
    q_melt <- melt(clustering_res, id.vars = "Sample")

    p <- ggplot(q_melt, aes(x = Sample, y = value, fill = variable)) +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
        labs(
            title = paste0("Population Structure (K=", best_k, ")"),
            y = "Ancestry coefficient", fill = "Cluster"
        )

    ggsave(file.path(analysis_results_dir, paste0(species_name, "_structure_plot.pdf")), p, width = 12, height = 6)
    cat("Structure plot saved.\n")
} else {
    cat("No structure detected (best K=1)\n")
}

cat("LEA analysis pipeline completed.\n")

# Load required libraries
if (!require("readr")) install.packages("readr", repos = "https://cran.r-project.org/")
if (!require("dplyr")) install.packages("dplyr", repos = "https://cran.r-project.org/")
if (!require("tidyr")) install.packages("tidyr", repos = "https://cran.r-project.org/")
if (!require("ggplot2")) install.packages("ggplot2", repos = "https://cran.r-project.org/")
if (!require("reshape2")) install.packages("reshape2", repos = "https://cran.r-project.org/")

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)

## starting with stacks catalogue generated against mbidens genome

## These loci have been generated using the algorithm of Maruki & Lynch 2015 that accounts for sequencing quality,
## allele balance and population-level allele frequencies. The genotyping process required the SNP and genotyping
## calling likelihood alpha value to be 0.01 (see 1.1.1), ensuring only high quality variants were called. This process
## takes these putative high quality genotype calls and filters them using the principles of Puritz et al (2018)


## this process uses the vcf.populations.snps.vcf file that has one row per SNP, but can have multiple lines per RAD locus.

# Step 1. Filter genotypes with genotype quality < 20 and minimun depth of 5
# Get VCF file path from environment variable or use fallback paths
vcf_file <- if (file.exists("populations.snps.vcf")) {
  "populations.snps.vcf"
} else if (!is.null(Sys.getenv("STACKS_DIR", unset = NA)) && !is.na(Sys.getenv("STACKS_DIR", unset = NA))) {
  file.path(Sys.getenv("STACKS_DIR"), "populations.snps.vcf")
} else if (file.exists("../results/stacks/pops/populations.snps.vcf")) {
  "../results/stacks/pops/populations.snps.vcf"
} else {
  # Full absolute path as fallback
  "/nesi/nobackup/uoa02626/Zyphius_demu/Analyse/BaseTest/bash-pipeline-project/results/stacks/pops/populations.snps.vcf"
}

# Verify the file exists before proceeding
if (!file.exists(vcf_file)) {
  stop(paste("VCF file not found at:", vcf_file))
}

cat("Using VCF file:", vcf_file, "\n")

# Create results directory for filtering outputs in the stacks results folder
if (!is.null(Sys.getenv("STACKS_DIR", unset = NA)) && !is.na(Sys.getenv("STACKS_DIR", unset = NA))) {
  # Use environment variable if available
  stacks_base_dir <- dirname(Sys.getenv("STACKS_DIR"))
  filter_results_dir <- file.path(stacks_base_dir, "filtered")
} else {
  # Fallback to relative path
  filter_results_dir <- "../results/stacks/filtered"
}

if (!dir.exists(filter_results_dir)) {
  dir.create(filter_results_dir, recursive = TRUE)
  cat("Created results directory:", filter_results_dir, "\n")
}

# Change to results directory for all outputs
setwd(filter_results_dir)
cat("Working directory set to:", getwd(), "\n")

system(paste("vcftools --vcf", vcf_file, "--out out.1 --minDP 5 --minGQ 20 --recode --recode-INFO-all"))

## count number of remaining loci
system("egrep -v '^#' out.1.recode.vcf > bald.out.1.recode.vcf")
## grab the third column, ie locus ID
system("cut bald.out.1.recode.vcf -f 3 > 1.loci.number")
## cut it so that you only have locus, not variable site, sort, list and count unique values
system("cut 1.loci.number -d ':' -f 1 | sort -u | wc -l")

# Step 2. Filter out sites that were made monomorphic by the previous filter
system("vcftools --vcf out.1.recode.vcf --maf 0.001 --out out.2 --recode --recode-INFO-all")

## count number of remaining loci
system("egrep -v '^#' out.2.recode.vcf > bald.out.2.recode.vcf")
## grab the third column, ie locus ID
system("cut bald.out.2.recode.vcf -f 3 > 2.loci.number")
## cut it so that you only have locus, not variable site, sort, list and count unique values
system("cut 2.loci.number -d ':' -f 1 | sort -u | wc -l")


# Step 3. Remove sites with more than 50% missing data
system("vcftools --vcf out.2.recode.vcf --out out.3 --max-missing 0.5 --recode --recode-INFO-all")

## count number of remaining loci
system("egrep -v '^#' out.3.recode.vcf > bald.out.3.recode.vcf")
## grab the third column, ie locus ID
system("cut bald.out.3.recode.vcf -f 3 > 3.loci.number")
## cut it so that you only have locus, not variable site, sort, list and count unique values
system("cut 3.loci.number -d ':' -f 1 | sort -u | wc -l")


# Cannot remove loci with extreme allele balance or filter by SNP quality as CH did
# but this is accounted for by the genotype calling algorithm

# Produce a file with missingness per individual
system("vcftools --vcf out.3.recode.vcf --out out.3 --missing-indv")

out_3_imiss <- read_tsv("out.3.imiss")

# Automatically detect the number of samples at this stage
n_samples_step3 <- nrow(out_3_imiss)
cat("Number of samples at step 3:", n_samples_step3, "\n")

species <- rep("Zca", n_samples_step3) # Use actual number of samples
missingness <- cbind(out_3_imiss, species)
missplot <- melt(missingness[, c(1, 5, 6)])
Species <- c("Zca")
missplot <- within(missplot, species <- factor(species, levels = Species))

b <- ggplot(data = missplot, aes(x = INDV, y = value, fill = species))
b <- b + geom_bar(stat = "identity")
b <- b + theme(
  axis.text.x = element_text(size = 6, angle = 90),
  axis.text.y = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_blank()
) +
  ggtitle("Missing data per sample") +
  labs(y = "Initial missing proportion") +
  ylim(0, 1)

pdf("missingness_per_sample.initial.pdf")
print(b)
dev.off()

# Select individuals with more than 50% missing data

miss_50 <- filter(out_3_imiss, F_MISS > 0.5) %>%
  select(INDV)

# Write the individuals to remove to a file
write_delim(miss_50, "remove.3.inds", col_names = FALSE)

# Step4: Remove individuals with >50% missing data
system("vcftools --vcf out.3.recode.vcf --out out.4 --remove remove.3.inds --recode --recode-INFO-all")

# Calculate site depth
system("vcftools --vcf out.4.recode.vcf --site-depth --out out.5")

site_depth_5 <- read_tsv("out.5.ldepth") %>%
  mutate(MEAN_DEPTH = SUM_DEPTH / 20) # 20 individus

# Step 5: Filter out loci with a mean site depth > 3x the overall mean
mean_site_depth_5 <- mean(site_depth_5$MEAN_DEPTH)
to_keep_5 <- filter(site_depth_5, MEAN_DEPTH < 3 * mean_site_depth_5)
mean_site_depth_5_filt <- mean(to_keep_5$MEAN_DEPTH)

# Plot the distribution again
pdf("meandepth_per_locus.S5initial.pdf")
qplot(to_keep_5$MEAN_DEPTH, main = "Mean Depth per locus (filtered)", xlab = "Mean Depth")
dev.off()

# Make a list of the sites to filter
to_filter_5 <- filter(site_depth_5, MEAN_DEPTH >= 3 * mean_site_depth_5) %>%
  select(CHROM, POS)

# Write the sites to remove to a file
write_delim(to_filter_5, "remove.5.sites", col_names = FALSE)

# Remove the sites with VCFtools
system("vcftools --vcf out.4.recode.vcf --out out.5 --exclude-positions remove.5.sites --recode --recode-INFO-all")

## count number of remaining loci
system("egrep -v '^#' out.5.recode.vcf > bald.out.5.recode.vcf")
## grab the third column, ie locus ID
system("cut bald.out.5.recode.vcf -f 3 > 5.loci.number")
## cut it so that you only have locus, not variable site, sort, list and count unique values
system("cut 5.loci.number -d ':' -f 1 | sort -u | wc -l")

## Step 6: create whitelist of loci

## at various steps sites have been removed - but also need to make sure the loci
## these sites are in have been removed. sometimes another site in a locus makes it through QC

## to do this we need a list of the chrom pos and stacks ID for all markers
## this is in populations.snps.vcf
system(paste("egrep -v '^#'", vcf_file, "> bald.out.all.recode.vcf"))
system("cut bald.out.all.recode.vcf -f 1,2,3 > loc.coord")

## read into R the chros, pos and ID
loc.pos <- read.delim("loc.coord", head = FALSE)
stacks.pos.all <- as.data.frame(loc.pos[, 3])
stacks.pos.all <- separate(stacks.pos.all, col = 1, sep = ":", into = c("s.loc", "s.pos"))
stacks.pos.allele.all <- cbind(loc.pos[, 1:2], stacks.pos.all)
colnames(stacks.pos.allele.all) <- c("CHROM", "POS", "s.loc", "s.pos")

loci.all <- unite(stacks.pos.allele.all, "chrom_pos", 1:2, sep = "_", remove = FALSE)

## now we get the info fo r
system("cut bald.out.5.recode.vcf -f 1,2,3 > loc.coord.initialWL")
loc.pos.initialWL <- read.delim("loc.coord.initialWL", head = FALSE)
stacks.pos.qc1 <- as.data.frame(loc.pos.initialWL[, 3])
stacks.pos.qc1 <- separate(stacks.pos.qc1, col = 1, sep = ":", into = c("s.loc", "s.pos"))
stacks.pos.qc1.1 <- cbind(loc.pos.initialWL[, 1:2], stacks.pos.qc1)
colnames(stacks.pos.qc1.1) <- c("CHROM", "POS", "s.loc", "s.pos")

loci.QC1 <- unite(stacks.pos.qc1.1, "chrom_pos", 1:2, sep = "_", remove = FALSE)

## want to know which loci in loci.all are NOT in loci.QC1

Dif <- setdiff(loci.all$chrom_pos, loci.QC1$chrom_pos)

# and get their stacks locus number
Dif.stacks <- loci.all[which(loci.all$chrom_pos %in% Dif), ]
Dif.stacks.loci <- unique(sort(Dif.stacks$s.loc))

## compare that with the initial good list
loci.QC1.stacks.loci <- unique(sort(loci.QC1$s.loc))

loci.QC2 <- loci.QC1[-which(Dif.stacks.loci %in% loci.QC1.stacks.loci), ]

length(unique(loci.QC2$s.loc))
whitelist <- as.data.frame(unique(loci.QC2$s.loc))

# generate white list
write_delim(unique(whitelist), "TruesGlobalQC.wl", col_names = FALSE)

## double check missingness before progressing - remove samples with >25% missingness
system("vcftools --vcf out.5.recode.vcf --out out.6 --missing-indv")

out_6_imiss <- read_tsv("out.6.imiss")

## Step 7
## re-look at sample-missingness with this whitelist

system("vcftools --vcf out.5.recode.vcf --out out.6 --missing-indv")

out_6_imiss <- read_tsv("out.6.imiss")

# Automatically detect the number of remaining samples
n_samples_final <- nrow(out_6_imiss)
cat("Number of samples after filtering:", n_samples_final, "\n")

species <- rep("Zca", n_samples_final) # Use actual number of remaining samples
missingness <- cbind(out_6_imiss, species)
missplot <- melt(missingness[, c(1, 5, 6)])
Species <- c("Zca")
missplot <- within(missplot, species <- factor(species, levels = Species))

c <- ggplot(data = missplot, aes(x = INDV, y = value, fill = species))
c <- c + geom_bar(stat = "identity")
c <- c + theme(
  axis.text.x = element_text(size = 6, angle = 90),
  axis.text.y = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.title = element_blank()
) +
  ggtitle("Missing data per sample") +
  labs(y = "Initial missing proportion") +
  ylim(0, 1)

print(c)

pdf("missingness_per_sample.step7.pdf")
print(c)
dev.off()



## get basic stats
system(paste("vcftools --vcf", vcf_file, "--out out.1 --missing-indv"))

out_1_imiss <- read_tsv("out.1.imiss")

## 28738 sites/loci

## but still some high coverage loci, so removed this and re-wrote the whitelist
system(paste("vcftools --vcf", vcf_file, "--site-depth --out out.2"))

site_depth_2 <- read_tsv("out.2.ldepth") %>%
  mutate(MEAN_DEPTH = SUM_DEPTH / 20) # 20 individus

mean_site_depth_2 <- mean(site_depth_2$MEAN_DEPTH)
to_keep_2 <- filter(site_depth_2, MEAN_DEPTH < 3 * mean_site_depth_2)
mean_site_depth_2_filt <- mean(to_keep_2$MEAN_DEPTH)

## filtered loci
# > mean(to_keep_2$MEAN_DEPTH)
# [1] 67.90641
# > sd(to_keep_2$MEAN_DEPTH)
# [1] 48.8953
# > range(to_keep_2$MEAN_DEPTH)
# [1]   1.804348 207.434783

# Make a list of the sites to filter
to_filter_2 <- filter(site_depth_2, MEAN_DEPTH >= 3 * mean_site_depth_2) %>%
  select(CHROM, POS)

# Write the sites to remove to a file
write_delim(to_filter_2, "remove.2.sites", col_names = FALSE)

# Remove the sites with VCFtools
system(paste("vcftools --vcf", vcf_file, "--out out.3 --exclude-positions remove.2.sites --recode --recode-INFO-all"))

## out.3 is what we want to work with next: identify discordant markers

#### extract whitelist from remaining loci
system("egrep -v '^#' out.4.recode.vcf > bald.out.4.recode.vcf")
system("cut bald.out.4.recode.vcf -f 3 > 4.loci.number")
system("cut 4.loci.number -d ':' -f 1 | sort -u > whitelistv2.loci")

#####################################################################################
## now we rerun populations, keeping only one sample in the replicate and using new whitelist
## output from here went to /nesi/nobackup/uoa02626/bw_ddrad/mbid_gen/pops_ssnp_whitelist
#####################################################################################

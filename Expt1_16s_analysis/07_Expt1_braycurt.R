library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)
library(ggtext)

# Data Preparation Section --------------------------------------------------

#set working directory
setwd("/Users/alyssadaigle/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s")

# Step 1: Prepare OTU matrix
otu_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE)
otu_mat <- as.matrix(otu_mat)

# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat <- as.matrix(tax_mat)

# Step 3: Prepare sample data  
samples_df <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_sample_metadata.tsv", header = TRUE)

# Data Filtering Section ----------------------------------------------------

# Subset OTU matrix up to "24-d6"
col_index <- which(colnames(otu_mat) == "24-d6")
otu_mat_subset <- otu_mat[, 1:col_index]

# Remove rows where all values are 0
otu_mat_subset_filtered <- otu_mat_subset[rowSums(otu_mat_subset) > 0, ]

# Subset the taxonomy matrix
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat_subset_filtered), ]

# Filter out specific taxa
tax_mat_filtered <- tax_mat_filtered[!(
    (tax_mat_filtered[, "Kingdom"] == "Bacteria" & 
         tax_mat_filtered[, "Phylum"] == "Cyanobacteriota" & 
         tax_mat_filtered[, "Class"] == "Chloroplast" & 
         tax_mat_filtered[, "Order"] == "Chloroplast" & 
         is.na(tax_mat_filtered[, "Family"])) | 
        is.na(tax_mat_filtered[, "Family"]) |
        tax_mat_filtered[, "Class"] == "Mitochondria"
), ]

# Ensure unique taxa names
rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))

# Filter OTU matrix based on taxonomy
otu_mat_final <- otu_mat_subset_filtered[rownames(otu_mat_subset_filtered) %in% rownames(tax_mat_filtered), ]

# Filter sample data for Expt 1 specific analysis
row_index <- which(rownames(samples_df) == "24-d6")
samples_df_filtered <- samples_df[1:row_index, , drop = FALSE]

# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
otu_mat_relative <- sweep(otu_mat_final, 2, colSums(otu_mat_final), `/`) * 100

# Filter taxa based on mean relative abundance (keep only those with >= 1%)
keep_taxa <- rownames(otu_mat_relative)[apply(otu_mat_relative, 1, function(x) any(x >= 1))]
otu_mat_filtered_counts <- otu_mat_final[keep_taxa, ]

# Filter taxonomy again based on newly calculated rel abundance in the OTU matrix
tax_mat_filtered_counts <- tax_mat_filtered[rownames(tax_mat_filtered) %in% rownames(otu_mat_filtered_counts), ]

# Phyloseq Object Creation --------------------------------------------------

# Step 12: Create phyloseq object
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(OTU, TAX, samples)

# bray curtis dissimilarity
# --------------------------------------------------

# Step 1: Extract and update metadata
meta <- as(sample_data(physeq), "data.frame") %>%
    rownames_to_column("SampleID") %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE) %>%
    column_to_rownames("SampleID")

# Step 2: Assign updated metadata back to the phyloseq object
sample_data(physeq) <- sample_data(meta)

# Step 3: Log-transform the counts (log(x + 1))
physeq_log <- transform_sample_counts(physeq, function(x) log1p(x))

# Step 4: Compute Bray-Curtis dissimilarity
bray_dist <- distance(physeq_log, method = "bray")

# Step 5: Perform PCoA
ordination <- ordinate(physeq_log, method = "PCoA", distance = bray_dist)

# Step 6: Plot ordination
plot_ordination(physeq_log, ordination, color = "micro", shape = "cyano") +
    geom_point(size = 3) +
    theme_cowplot()



# MCMCglmm to see if treatment predicts the axes
# --------------------------------------------------

# Step 1: Extract the PCoA scores
pcoa_scores <- as.data.frame(ordination$vectors)

# Step 2: Add sample metadata
pcoa_scores$SampleID <- rownames(pcoa_scores)

pcoa_scores <- merge(pcoa_scores, samples_df_filtered, by.x = "SampleID", by.y = "row.names", all.x = TRUE)

# Step 3: separate treatment into components
pcoa_scores <- pcoa_scores %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)


mod1 <- MCMCglmm(Axis.2 ~ cyano + micro, data = pcoa_scores, nitt = 11000, burnin = 1000, thin = 10, verbose = FALSE)
summary(mod1)

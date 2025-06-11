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

# PCA Analysis Section -------------------------------------------------------
otu_data <- otu_table(physeq)
otu_data_log <- log1p(otu_data)

# Perform PCA using prcomp()
pca_results <- prcomp(t(otu_data_log), scale. = TRUE)
eig_vals <- eigenvals(pca_results)  # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100  # Variance explained

# Extract PCA scores for the first three components (PC1, PC2, and PC3)
pca_scores <- as.data.frame(pca_results$x[, 1:3])

# If you want to merge with your sample metadata:
pca_scores$SampleID <- rownames(pca_scores)

pca_scores <- merge(pca_scores, samples_df_filtered, by.x = "SampleID", by.y = "row.names", all.x = TRUE)

# separate treatment into components
pca_scores <- pca_scores %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)

# assign colors to microbes
custom_colors <- c("#AA4499", "#DDCC77", "#88CCEE", "#117733")

# Step 4: Visualize PCA
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = micro, shape = cyano)) +
    geom_point(size = 2.5) +
    scale_color_manual(name = "Microbiome Source",
                       values = custom_colors,
                       breaks = c("H", "KF", "ODR", "N"),
                       labels = c("Home", "Kingman Farm", "Dairy Farm", "Uninoculated")) +
    scale_shape_manual(name = "*M. aeruginosa* Spike",
                       values = c("Y" = 16, "N" = 17),
                       labels = c("Y" = "Yes", "N" = "No")) +
    labs(x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
         y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = "")) +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme_cowplot() +
    theme(legend.title = element_markdown(size = 10),
          legend.text = element_text(size = 8), 
          strip.text = element_text(size = 5.7))


pca_plot 

#save plot
ggsave("Expt1_pca_plot.jpg", pca_plot, width = 6, height = 4)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_pca_plot.jpg", pca_plot, width = 6, height = 4)



# Calculate mean and standard deviation for each group of microbiome and cyano
pca_stats <- pca_scores %>%
    group_by(micro, cyano) %>%
    summarise(mean_PC1 = mean(PC1),
              sd_PC1 = sd(PC1),
              .groups = 'drop')

dodge_width <- 1

# Create the plot
pc1_micro_plot <- ggplot(pca_stats, aes(x = factor(micro, levels = c("N", "H", "KF", "ODR")), y = mean_PC1, color = cyano)) +
    geom_jitter(data = pca_scores,
                aes(x = factor(micro, levels = c("N", "H", "KF", "ODR")), 
                    y = PC1, color = cyano),
                size = 1.5,
                shape = 1,
                position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = dodge_width),
                inherit.aes = FALSE) +
    geom_point(size = 2, shape = 19, position = position_dodge(width = dodge_width)) + 
    geom_errorbar(aes(ymin = mean_PC1 - sd_PC1, ymax = mean_PC1 + sd_PC1),
                  width = 0.2,
                  position = position_dodge(width = dodge_width)) +
    scale_color_manual(name = "*M. aeruginosa* Spike",
                       values = c("Y" = "aquamarine4", "N" = "black"),
                       labels = c("Y" = "Yes", "N" = "No")) +
    scale_x_discrete(labels = c("N" = "Uninoculated", 
                                "H" = "Home", 
                                "KF" = "Kingman Farm", 
                                "ODR" = "Dairy Farm")) +
    labs(x = "Microbiome Source",
         y = paste("PC1 (", round(varexpl[1], 2), "%)", sep = "")) +
    theme_cowplot() +
    theme(legend.title = element_markdown(size = 8),
          legend.text = element_text(size = 7),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 9))

# Print the plot
print(pc1_micro_plot)


ggsave("Expt1_pc1_micro_plot.jpg", pc1_micro_plot, width = 5, height = 3)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_pc1_micro_plot.jpg", pc1_micro_plot, width = 5, height = 3)


# which taxa are driving the PCA

# Load species (OTU) scores
otu_loadings <- scores(pca_results, display = "species") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("OTU")

# Extract taxonomy table
taxonomy_df <- tax_table(physeq) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("OTU") %>%
    select(OTU, Family)

# Function to get top N OTUs for a given PC
get_top_otus <- function(df, pc_col, n = 10) {
    df %>%
        arrange(desc(abs(.data[[pc_col]]))) %>%
        slice_head(n = n) %>%
        left_join(taxonomy_df, by = "OTU") %>%
        select(OTU, !!pc_col, Family)
}

# Get top 10 for PC1 and PC2
top10_PC1_df <- get_top_otus(otu_loadings, "PC1")
top10_PC2_df <- get_top_otus(otu_loadings, "PC2")


#  linear model to see that PC1 is significantly affected by treatment

pca_scores$treatment <- as.factor(pca_scores$treatment)

lm2_PC1 <- MCMCglmm(PC1 ~ -1 + cyano + micro + geno, data = pca_scores, nitt = 101000, burnin = 1000, thin = 10, verbose = FALSE)
summary(lm2_PC1)

lm_PC2 <- MCMCglmm(PC2 ~ -1 + cyano + micro + geno, data = pca_scores, nitt = 101000, burnin = 1000, thin = 10, verbose = FALSE)
summary(lm_PC2)

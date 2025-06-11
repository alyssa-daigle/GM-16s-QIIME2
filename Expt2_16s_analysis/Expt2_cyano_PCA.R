library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)
library(viridis)
library(pals)
library(Polychrome)
library(rcartocolor)
library(ggtext)
library(patchwork)

# Data Preparation Section --------------------------------------------------

#set working directory
setwd("/Users/alyssadaigle/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s")

# Step 1: Prepare OTU matrix
otu_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_otu_matrix.tsv", header = TRUE, check.names = FALSE)
otu_mat <- as.matrix(otu_mat)

# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat <- as.matrix(tax_mat)

# Step 3: Prepare sample data  
samples_df <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_sample_metadata.tsv", header = TRUE)


# Data Filtering Section ----------------------------------------------------

# Subset OTU matrix after MP-1-1
col_index <- which(colnames(otu_mat) == "MP-1-1")
otu_mat_subset <- otu_mat[, col_index:ncol(otu_mat)]
otu_mat_subset <- otu_mat_subset[, !grepl("Expt1|soil|blank", colnames(otu_mat_subset))]
otu_mat_subset <- otu_mat_subset[rowSums(otu_mat_subset) > 0, ]

# Filter taxonomy to match OTU
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat_subset), ]

# Filter for Cyanophyceae 
tax_mat_cyano <- tax_mat_filtered %>%
    as.data.frame() %>%
    filter(Class == "Cyanophyceae") %>%
    drop_na(Family) %>%
    filter(Family != c("Unclassified", "Unassigned"))

# Keep only rows in OTU table that match Cyanophyceae
otu_mat_cyano <- otu_mat_subset[rownames(otu_mat_subset) %in% rownames(tax_mat_cyano), ]

# Convert to long format and merge taxonomy early
otu_long_cyano <- as.data.frame(otu_mat_cyano) %>%
    rownames_to_column("FeatureID") %>%
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") %>%
    left_join(tax_mat_cyano %>% rownames_to_column("FeatureID"), by = "FeatureID")

# Add metadata (filtered only once)
samples_df_filtered <- samples_df %>%
    filter(!grepl("Expt1|soil|blank", rownames(samples_df))) %>%
    mutate(SampleID = rownames(.))

# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
otu_mat_relative <- sweep(otu_mat_cyano, 2, colSums(otu_mat_cyano), `/`) * 100

# Filter taxa based on mean relative abundance
keep_taxa <- rownames(otu_mat_relative)[apply(otu_mat_relative, 1, function(x) any(x >= 1))]
otu_mat_filtered_counts <- otu_mat_cyano[keep_taxa, ]

# Filter taxonomy again based on newly calculated rel abundance in the OTU matrix
tax_mat_filtered_counts <- tax_mat_cyano[rownames(tax_mat_cyano) %in% rownames(otu_mat_filtered_counts), ]

# Phyloseq Object Creation --------------------------------------------------

# Step 12: Create phyloseq object
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
tax_mat_filtered_counts <- as.matrix(tax_mat_filtered_counts)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(OTU, TAX, samples)

# PCA Analysis ------------------------------------------------

otu_data <- otu_table(physeq)
otu_data_log <- log1p(otu_data)

# Step 2: Perform PCA
pca_results <- rda(t(otu_data_log))  # Transpose the data for PCA
eig_vals <- eigenvals(pca_results)  # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100  # Variance explained

# Step 3: Extract PCA scores
pca_scores <- as.data.frame(scores(pca_results, display = "sites"))
pca_scores$SampleID <- rownames(pca_scores)
pca_scores <- merge(pca_scores, samples_df_filtered, by.x = "SampleID", by.y = "row.names", all.x = TRUE)

# Extract pond category (first part before second dash)
pca_scores$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", pca_scores$treatment)

# Categorize treatment as 'Water' or 'Duckweed'
pca_scores$treatment_category <- ifelse(grepl("DW", pca_scores$treatment), "Duckweed", "Water")

# Mapping pond names to full names
pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson Farm\nPond 1",
    "TF-2" = "Thompson Farm\nPond 2",
    "UM-1" = "Upper Mill Pond")

# Recode pond names using pond_name_mapping
pca_scores$pond_full_name <- recode(pca_scores$pond, !!!pond_name_mapping)


# Define the custom color palette manually
custom_colors <- c("#44AA99", "#CC6677", "#88CCEE", "#117733", "#DDCC77", "#332288")

# Apply the selected colors in ggplot
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pond_full_name, shape = treatment_category)) +
    geom_point(size = 3.5) +
    scale_color_manual(values = custom_colors) +  # Use manually selected colors
    labs(x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
         y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = ""),
         color = "Green Manure Source", 
         shape = "Sample Type") +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2)) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 5.7))

print(pca_plot)


#save plot
ggsave("expt2_cyano_pca_plot.jpg", pca_plot, width = 6, height = 4)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/expt2_cyano_pca_plot.jpg", pca_plot, width = 6, height = 4)

# Duckweed-specific PCA plot --------------------------

# filter the physeq object to just duckweed samples and perform new PCA
physeq_dw <- subset_samples(physeq, grepl("DW", treatment))
# Perform PCA on duckweed samples
dw_otu_data <- otu_table(physeq_dw)
dw_otu_data_log <- log1p(dw_otu_data)

# Step 2: Perform PCA
pca_results <- rda(t(dw_otu_data_log))  # Transpose the data for PCA
eig_vals <- eigenvals(pca_results)  # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100  # Variance explained

# Step 3: Extract PCA scores
pca_scores <- as.data.frame(scores(pca_results, display = "sites"))
pca_scores$SampleID <- rownames(pca_scores)
pca_scores <- merge(pca_scores, samples_df_filtered, by.x = "SampleID", by.y = "row.names", all.x = TRUE)

# Mapping pond names to full names
pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson Farm\nPond 1",
    "TF-2" = "Thompson Farm\nPond 2",
    "UM-1" = "Upper Mill Pond")

#separate treatment into pond
pca_scores$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", pca_scores$treatment)

# Define the custom color palette manually
custom_colors <- c("#44AA99", "#CC6677", "#88CCEE", "#117733", "#DDCC77", "#332288")

pca_scores$pond <- recode(pca_scores$pond, !!!pond_name_mapping)

# Apply the selected colors in ggplot
dw_pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pond)) +
    geom_point(size = 3.5) +
    scale_color_manual(values = custom_colors, labels = pond_name_mapping) +  # Use manually selected colors
    labs(x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
         y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = ""),
         color = "Green Manure Source") +
    theme_cowplot() +
    theme(
        legend.position = "none",
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 5.7))

# Combine plots, keeping the legend only from the first
combined_plot <- pca_plot + dw_pca_plot + plot_layout(ncol = 1, guides = "collect") & theme(legend.position = "right")

# Display the combined plot
print(combined_plot)

ggsave("expt2_DWcyano_pca_plot.jpg", combined_plot, width = 6, height = 8)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/expt2_DWcyano_pca_plot.jpg", combined_plot, width = 6, height = 6)

# -----------------------------------------------------------
# which taxa are driving the PCA
# also linear model to see that PC1 is significantly affected by treatmemnt

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

# models -----------------------------------

# Ensure treatment is a factor (if not already)
pca_scores$treatment_category <- as.factor(pca_scores$treatment_category)

PC1 <- MCMCglmm(PC1 ~ -1 + treatment_category:pond, data = pca_scores, nitt = 101000, burnin = 1000, thin = 10, verbose = FALSE)
summary(PC1)

PC2 <- MCMCglmm(PC2 ~ treatment_category:pond, data = pca_scores, nitt = 13000, burnin = 3000, thin = 10, verbose = FALSE)
summary(PC2)



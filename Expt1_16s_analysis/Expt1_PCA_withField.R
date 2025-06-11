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

# Step 1: Prepare OTU matrix
otu_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_otu_matrix.tsv", header = TRUE, row.names = 1, check.names = FALSE)
otu_mat <- as.matrix(otu_mat)

# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat <- as.matrix(tax_mat)

# Step 3: Prepare sample data  
samples_df <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/processed_sample_metadata.tsv", header = TRUE)

# Data Filtering Section ----------------------------------------------------

# Find the index of the "24-d6" column
col_index <- which(colnames(otu_mat) == "24-d6")

# Get the base subset (columns from 1 to col_index)
base_cols <- 1:col_index

# Get additional columns that end in "-Expt1", but exclude "W-soil-Expt1"
extra_cols <- which(grepl("-Expt1$", colnames(otu_mat)) & colnames(otu_mat) != "W-soil-Expt1")

# Combine the column indices and ensure uniqueness
final_cols <- sort(unique(c(base_cols, extra_cols)))

# Subset the matrix
otu_mat_subset <- otu_mat[, final_cols]

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

# Find the index of the "24-d6" row
row_index <- which(rownames(samples_df) == "24-d6")

# Get base rows: 1 through "24-d6"
base_rows <- 1:row_index

# Get additional rows that end in "-Expt1", excluding "W-soil-Expt1"
extra_rows <- which(grepl("-Expt1$", rownames(samples_df)) & rownames(samples_df) != "W-soil-Expt1")

# Combine and sort unique row indices
final_rows <- sort(unique(c(base_rows, extra_rows)))

# Subset the dataframe
samples_df_filtered <- samples_df[final_rows, , drop = FALSE]

# Find the index for "24-d6"
cutoff_row <- which(rownames(samples_df_filtered) == "24-d6")

# Separate the data into two parts
df_tidy <- samples_df_filtered[1:cutoff_row, ]
df_rest <- samples_df_filtered[(cutoff_row + 1):nrow(samples_df_filtered), ]

# Separate the treatment column in the top part only
df_tidy <- df_tidy %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)

# Combine the two parts back together
samples_df_filtered <- bind_rows(df_tidy, df_rest)

samples_df_filtered <- samples_df_filtered %>%
    mutate(inoculum = case_when(
        micro == "N" ~ "Uninoculated",
        micro == "H" ~ "Home",
        micro == "ODR" ~ "Dairy Farm",
        micro == "KF" ~ "Kingman Farm",
        TRUE ~ NA_character_
    )) |> 
    select(-geno, -cyano, -micro)

samples_df_filtered <- samples_df_filtered %>%
    mutate(inoculum = case_when(
        grepl("W-Expt1", treatment)   ~ "Field",
        grepl("Dr-Expt1", treatment)  ~ "Field",
        grepl("LR-Expt1", treatment)  ~ "Field",
        grepl("UM-Expt1", treatment)  ~ "Field",
        grepl("MP-Expt1", treatment)  ~ "Field",
        grepl("TF-Expt1", treatment)  ~ "Field",
        grepl("KF-Expt1", treatment)  ~ "Field",
        grepl("ODR-Expt1", treatment) ~ "Field",
        TRUE ~ inoculum  # keep existing values for all other treatments
    ))

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
# Separate treatment only in first 192 rows
first_part <- pca_scores[1:192, ] %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)

# For rows after 192, add the same columns but fill with NA
second_part <- pca_scores[193:nrow(pca_scores), ]
second_part$geno <- NA
second_part$cyano <- "N"
second_part$micro <- NA

# Recombine
pca_scores <- bind_rows(first_part, second_part)

# assign colors to microbes
custom_colors <- c("#AA4499", "#DDCC77", "#88CCEE", "#117733", "orange")

#set inoculum levels 
pca_scores$inoculum <- factor(pca_scores$inoculum, levels = c("Home", "Kingman Farm", "Dairy Farm", "Uninoculated", "Field"))

# Step 4: Visualize PCA
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = inoculum, shape = cyano)) +
    geom_point(size = 2.5) +
    scale_color_manual(name = "Microbiome Source",
                       values = custom_colors)+
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
                                              
lm_PC1 <- MCMCglmm(PC1 ~ -1 + inoculum, data = pca_scores, nitt = 11000, burnin = 1000, thin = 10, verbose = FALSE)
summary(lm_PC1)

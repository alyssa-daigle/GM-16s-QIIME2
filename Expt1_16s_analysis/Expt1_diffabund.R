library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)
library(ggtext)
library(corncob)

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

samples_df_filtered$treatment <- as.character(samples_df_filtered$treatment)

# Separate treatment into geno, cyano, and micro columns
samples_df_filtered <- samples_df_filtered %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)
   # mutate(micro = if_else(micro == "H", paste0("H_", geno), micro))

# Make sure 'micro' is a factor with microN as the reference level
# samples_df_filtered$micro <- factor(samples_df_filtered$micro, levels = c("N", "H_DR", "H_LR", "H_M", "H_TF", "H_UM", "H_W", "KF", "ODR"))

samples_df_filtered <- samples_df_filtered %>%
    mutate(
        inoc_status = factor(
            if_else(cyano == "N", "uninoculated", "inoculated"),
            levels = c("uninoculated", "inoculated")
        )
    )

## Micro = N=specific analysis ---------------------------------------------
microN_df <- samples_df_filtered %>%
    filter(micro == "N")

# Create phyloseq components
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(microN_df)
microN_physeq <- phyloseq(OTU, TAX, samples)

# Agglomerate taxa at the Family level
microN_physeq_family <- tax_glom(microN_physeq, taxrank = "Family")

# Convert tax_table to data frame
tax_df <- as.data.frame(tax_table(microN_physeq_family))

# Identify taxa to remove
families_to_remove <- c("NS11-12 marine group", "Mitochondria", "Incertae Sedis")

# Identify taxa to keep (i.e., not in the families_to_remove)
taxa_to_keep <- rownames(tax_df)[!(tax_df$Family %in% families_to_remove)]

# Prune the phyloseq object
microN_physeq_filtered <- prune_taxa(taxa_to_keep, microN_physeq_family)

# Run Differential Abundance Testing with corncob 

diff_test <- differentialTest(
    formula = ~ inoc_status,
    phi.formula = ~ inoc_status,
    formula_null = ~ 1,
    phi.formula_null = ~ inoc_status,
    test = "Wald",
    boot = FALSE,
    data = microN_physeq_filtered,
    fdr_cutoff = 0.05,
    full_output = FALSE)

# Extract plotting data just for those taxa
plot_data <- plot(diff_test, level = "Family", data_only = TRUE)

# Step 2: Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data |> 
    drop_na() |>
    group_by(variable) |> 
    mutate(taxa = forcats::fct_reorder(taxa, x)) |> 
    ungroup() 

microN_plot <- ggplot(plot_data, aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(title = "Uninoculated",
         x = ""
    ) +
    theme_bw() +
    theme(
        title = element_text(size = 6, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, hjust = 0.5), 
        legend.position = "none",
        axis.title.y = element_blank())

microN_plot

## Micro = H-specific analysis ---------------------------------------------
microH_df <- samples_df_filtered %>%
    filter(micro == "H")

# Create phyloseq components
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(microH_df)
microH_physeq <- phyloseq(OTU, TAX, samples)

# Agglomerate taxa at the Family level
microH_physeq_family <- tax_glom(microH_physeq, taxrank = "Family")

# Convert tax_table to data frame
tax_df <- as.data.frame(tax_table(microH_physeq_family))

# Identify taxa to remove
families_to_remove <- c("NS11-12 marine group", "Mitochondria", "Incertae Sedis")

# Identify taxa to keep (i.e., not in the families_to_remove)
taxa_to_keep <- rownames(tax_df)[!(tax_df$Family %in% families_to_remove)]

# Prune the phyloseq object
microH_physeq_filtered <- prune_taxa(taxa_to_keep, microH_physeq_family)

# Run Differential Abundance Testing with corncob 
micrH_diff_test <- differentialTest(
    formula = ~ inoc_status,
    phi.formula = ~ inoc_status,
    formula_null = ~ 1,
    phi.formula_null = ~ inoc_status,
    test = "Wald",
    boot = FALSE,
    data = microH_physeq_filtered,
    fdr_cutoff = 0.05,
    full_output = FALSE)

# Extract plotting data just for those taxa
plot_data <- plot(micrH_diff_test, level = "Family", data_only = TRUE)

# Step 2: Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data |> 
    drop_na() |>
    group_by(variable) |> 
    mutate(taxa = forcats::fct_reorder(taxa, x)) |> 
    ungroup() 

microH_plot <- ggplot(plot_data, aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(title = "Home microbiome",
         x = ""
    ) +
    theme_bw() +
    theme(
        title = element_text(size = 6, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, hjust = 0.5), 
        legend.position = "none",
        axis.title.y = element_blank())
microH_plot

## Micro = KF-specific analysis ---------------------------------------------
microKF_df <- samples_df_filtered %>%
    filter(micro == "KF")

# Create phyloseq components
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(microKF_df)
microKF_physeq <- phyloseq(OTU, TAX, samples)

# Agglomerate taxa at the Family level
microKF_physeq_family <- tax_glom(microKF_physeq, taxrank = "Family")

# Convert tax_table to data frame
tax_df <- as.data.frame(tax_table(microKF_physeq_family))

# Identify taxa to remove
families_to_remove <- c("NS11-12 marine group", "Mitochondria", "Incertae Sedis")

# Identify taxa to keep (i.e., not in the families_to_remove)
taxa_to_keep <- rownames(tax_df)[!(tax_df$Family %in% families_to_remove)]

# Prune the phyloseq object
microKF_physeq_filtered <- prune_taxa(taxa_to_keep, microKF_physeq_family)

# Run Differential Abundance Testing with corncob 
microKF_diff_test <- differentialTest(
    formula = ~ inoc_status,
    phi.formula = ~ inoc_status,
    formula_null = ~ 1,
    phi.formula_null = ~ inoc_status,
    test = "Wald",
    boot = FALSE,
    data = microKF_physeq_filtered,
    fdr_cutoff = 0.05,
    full_output = FALSE)

# Extract plotting data just for those taxa
plot_data <- plot(microKF_diff_test, level = "Family", data_only = TRUE)

# Step 2: Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data |> 
    drop_na() |>
    group_by(variable) |> 
    mutate(taxa = forcats::fct_reorder(taxa, x)) |> 
    ungroup() 

microKF_plot <- ggplot(plot_data, aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(title = "Kingman Farm microbiome",
         x = ""
    ) +
    theme_bw() +
    theme(
        title = element_text(size = 6, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, hjust = 0.5), 
        legend.position = "none",
        axis.title.y = element_blank())
microKF_plot

## Micro = ODR-specific analysis --------------------------------------------
microODR_df <- samples_df_filtered %>%
    filter(micro == "ODR")

# Create phyloseq components
OTU <- otu_table(otu_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(microODR_df)
microODR_physeq <- phyloseq(OTU, TAX, samples)

# Agglomerate taxa at the Family level
microODR_physeq_family <- tax_glom(microODR_physeq, taxrank = "Family")

# Convert tax_table to data frame
tax_df <- as.data.frame(tax_table(microODR_physeq_family))

# Identify taxa to remove
families_to_remove <- c("NS11-12 marine group", "Mitochondria", "Incertae Sedis")

# Identify taxa to keep (i.e., not in the families_to_remove)
taxa_to_keep <- rownames(tax_df)[!(tax_df$Family %in% families_to_remove)]

# Prune the phyloseq object
microODR_physeq_filtered <- prune_taxa(taxa_to_keep, microODR_physeq_family)

# Run Differential Abundance Testing with corncob 
microODR_diff_test <- differentialTest(
    formula = ~ inoc_status,
    phi.formula = ~ inoc_status,
    formula_null = ~ 1,
    phi.formula_null = ~ inoc_status,
    test = "Wald",
    boot = FALSE,
    data = microODR_physeq_filtered,
    fdr_cutoff = 0.05,
    full_output = FALSE)

# Extract plotting data just for those taxa
plot_data <- plot(microODR_diff_test, level = "Family", data_only = TRUE)

# Step 2: Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data |> 
    drop_na() |>
    group_by(variable) |> 
    mutate(taxa = forcats::fct_reorder(taxa, x)) |> 
    ungroup() 

microODR_plot <- ggplot(plot_data, aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(title = "Dairy Farm microbiome",
         x = "",
         color = "Direction \nof Effect"
    ) +
    theme_bw() +
    theme(
        title = element_text(size = 6, face = "bold"),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, hjust = 0.5), 
        legend.position = "none", 
        axis.title.y = element_blank())

microODR_plot

## save all plots together
combined_plot <- plot_grid(microN_plot, microH_plot, microKF_plot, microODR_plot, ncol = 4)
combined_plot

ggsave("Expt1_diffabund_plot_combined.jpg", width = 10, height = 4, dpi = 500)
# ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_diffabund_plot_combined.jpg", width = 5, height = 6, dpi = 500)
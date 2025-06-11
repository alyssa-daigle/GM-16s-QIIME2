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

#  Remove rows where all values are 0
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

# Filter sample data for Expt 2 specific analysis
row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[!grepl("Expt1|soil|blank", rownames(samples_df_filtered)), , drop = FALSE]

# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
# fix this so its col sums
otu_mat_relative <- sweep(otu_mat_final, 2, colSums(otu_mat_final), `/`) * 100

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

# Remove taxa where Class is unassigned or Chloroplast
physeq <- subset_taxa(physeq, !(Class %in% c("unassigned", "Chloroplast")))

# Extract sample data
sample_data_df <- data.frame(phyloseq::sample_data(physeq))

# Create new column based on treatment values
sample_data_df$treatment_type <- factor(
    ifelse(grepl("DW$", sample_data_df$treatment), "duckweed", "water"),
    levels = c("water", "duckweed")
)

# Update the phyloseq object
phyloseq::sample_data(physeq) <- sample_data_df

physeq_class <- tax_glom(physeq, taxrank = "Class")

# Run the corncob test on the filtered phyloseq object
diff_test <- differentialTest(
    formula = ~ treatment_type,,
    phi.formula = ~ treatment_type,
    formula_null = ~ 1,
    phi.formula_null = ~ treatment_type,
    test = "Wald",
    boot = FALSE,
    data = physeq_class,
    fdr_cutoff = 0.05)

# Step 1: Extract the plot data
plot_data <- plot(diff_test, level = "Class", data_only = TRUE)
plot_data$p_fdr <- diff_test$p_fdr[plot_data$taxa]

# Step 2: Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data %>%
    group_by(variable) %>%
    mutate(taxa = forcats::fct_reorder(taxa, x)) %>%
    ungroup()

print(plot_data)

ggplot(plot_data, aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(
        x = "log(Odds Ratio) relative to Water Samples",
        y = "Bacterial Class",
        color = "Direction of Effect"
    ) +
    theme_bw() +
    theme(
        axis.text.y = element_text(size = 8, margin = margin(l = 0.4, r=0.1, unit = "cm")),
        axis.text.x = element_text(angle = 0, hjust = 0.5))


ggsave("Expt2_diffabund_plot.jpg", width = 5, height = 4, dpi = 500)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt2_diffabund_plot.jpg", width = 5, height =4, dpi = 500)



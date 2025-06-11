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
colnames(samples_df)[2] <- "SampleID"

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

# Shannon Diversity analysis ---------------------------------------------------
shannon_div <- diversity(otu_mat_final, index = "shannon", MARGIN = 2)

shannon_div_df <- data.frame(
    SampleID = colnames(otu_mat_final),  # Sample names are the column names
    ShannonDiversity = shannon_div
)

# Merge with sample metadata
shannon_div_df <- merge(shannon_div_df, samples_df_filtered, by = "SampleID")

# Step 1: Define and preprocess the variables
# Separate the treatment into geno, cyano, and micro
shannon_div_df <- shannon_div_df %>%
    separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE)

# Reorder the 'micro' factor to match the label order
shannon_div_df$micro <- factor(shannon_div_df$micro, levels = c("N", "H", "KF", "ODR"))

# Update the cyano labels in the data
shannon_div_df$cyano <- factor(shannon_div_df$cyano, levels = c("N", "Y"), labels = c("No", "Yes"))

# Create a new 'x_axis' based on the 'geno' and 'inoc' components
shannon_div_df$group <- sub("-[YN]-.*", "", shannon_div_df$treatment)  # Extract group (e.g., DR, LR, etc.)
shannon_div_df$inoc <- sub("^.*?-([YN])-.*", "\\1", shannon_div_df$treatment)  # Extract inoculation (Y or N)
shannon_div_df$x_axis <- paste(shannon_div_df$group, shannon_div_df$inoc, sep = "_")  # Combine group and inoc

# Step 2: Define levels for x_axis and labels for the x-axis
x_levels <- c("DR_Y", "DR_N", "LR_Y", "LR_N", "M_Y", "M_N", "TF_Y", "TF_N", "UM_Y", "UM_N", "W_Y", "W_N")
inoc_labels <- c("Y", "N", "Y", "N", "Y", "N", "Y", "N", "Y", "N", "Y", "N")

# Ensure the x_axis is a factor with the predefined levels
shannon_div_df$x_axis <- factor(shannon_div_df$x_axis, levels = x_levels)

micro_labels <- c("N" = "Uninoculated",
                  "H" = "Home microbiome",
                  "KF" = "Kingman Farm microbiome",
                  "ODR" = "Dairy Farm microbiome")

cyano_colors <- c("No" = "black", "Yes" = "aquamarine4")

# Set the labels for the x-axis, keeping Y/N only (no group labels like DR, LR, etc.)
x_labels <- setNames(inoc_labels, x_levels)

# Step 3: Plot
dodge_width <- 0.5  # Define dodge width

# Create the updated Shannon Diversity plot
Expt1_shan_div <- ggplot(shannon_div_df, aes(x = geno, y = ShannonDiversity, color = cyano)) +
    stat_summary(fun = "mean", geom = "point", size = 1.5, shape = 19, 
                 position = position_dodge(width = dodge_width)) +
    stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3, 
                 position = position_dodge(width = dodge_width)) + 
    theme_cowplot() +
    facet_wrap(~micro, scales = "free_x", ncol = 4, labeller = labeller(micro = micro_labels)) +
    theme(
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "right",
        legend.title = element_markdown(size = 8),
        legend.text = element_markdown(size =8)) +
    scale_color_manual(values = cyano_colors) +
    scale_x_discrete(labels = x_labels) +
    labs(x = "Duckweed Genotype", y = "Shannon Diversity Index", color = "*M. aeruginosa* spike")

Expt1_shan_div

# Save the plot
ggsave("Expt1_shannon_diversity_barplot.jpg", Expt1_shan_div, width = 12, height = 4)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_shannon_diversity_barplot.jpg", Expt1_shan_div, width = 12, height = 4)

# Linear model section ------------------------------------------------

#MCMCglmm 
mcmc_data <- shannon_div_df %>%
    select(ShannonDiversity, treatment, geno, cyano, micro) %>%
    mutate(geno = as.factor(geno),
           cyano = as.factor(cyano),
           micro = as.factor(micro))



mod1 <- MCMCglmm(ShannonDiversity ~  -1 + micro, data=mcmc_data,verbose=F, nitt = 101000, thin = 10, burnin = 1000) 
summary(mod1)

post_summ <- summary(mod1)$solutions
ci_df <- as.data.frame(post_summ)
ci_df$Effect <- rownames(ci_df)
names(ci_df)[c(1,2,3)] <- c("PostMean", "Lower95CI", "Upper95CI")

ggplot(ci_df, aes(x = Effect, y = PostMean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
    theme_minimal() +
    labs(title = "Posterior Means and 95% Credible Intervals",
         y = "Posterior Mean ± 95% CI",
         x = "Effect") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))


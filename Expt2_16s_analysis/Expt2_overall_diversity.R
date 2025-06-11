library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(tibble)
library(MCMCglmm)

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

# Step 4: Subset OTU matrix after MP-1-1
col_index <- which(colnames(otu_mat) == "MP-1-1")
otu_mat_subset <- otu_mat[, col_index:ncol(otu_mat)]
otu_mat_subset <- otu_mat_subset[, !grepl("Expt1|soil|blank", colnames(otu_mat_subset))]


# Step 5: Remove rows where all values are 0
otu_mat_subset_filtered <- otu_mat_subset[rowSums(otu_mat_subset) > 0, ]

# Step 6: Subset the taxonomy matrix
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat_subset_filtered), ]

# Step 7: Filter out specific taxa
tax_mat_filtered <- tax_mat_filtered[!(
    (tax_mat_filtered[, "Kingdom"] == "Bacteria" & 
         tax_mat_filtered[, "Phylum"] == "Cyanobacteriota" & 
         tax_mat_filtered[, "Class"] == "Chloroplast" & 
         tax_mat_filtered[, "Order"] == "Chloroplast" & 
         is.na(tax_mat_filtered[, "Family"])) | 
        is.na(tax_mat_filtered[, "Family"]) |
        tax_mat_filtered[, "Class"] == "Mitochondria"
), ]

# Step 8: Ensure unique taxa names
rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))

# Step 9: Filter OTU matrix based on taxonomy
otu_mat_final <- otu_mat_subset_filtered[rownames(otu_mat_subset_filtered) %in% rownames(tax_mat_filtered), ]

# Step 11: Filter sample data for Expt 2 specific analysis
row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[!grepl("Expt1|soil", rownames(samples_df_filtered)), , drop = FALSE]

otu_mat_subset <- otu_mat_subset[, !grepl("Expt1|soil|blank", colnames(otu_mat_subset))]

# Diversity Analysis Section ------------------------------------------------

# Step 1: Calculate alpha diversity (Shannon index) for each sample
shannon_div <- diversity(otu_mat_final, index = "shannon", MARGIN = 2)

shannon_div_df <- data.frame(
    SampleID = colnames(otu_mat_final),
    treatment = colnames(otu_mat_final),
    ShannonDiversity = shannon_div)

# Merge with sample metadata
shannon_div_df <- merge(shannon_div_df, samples_df_filtered, by = "treatment")

# extract treatment components
shannon_div_df$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", shannon_div_df$treatment)

# Categorize treatment as 'Water' or 'Duckweed'
shannon_div_df$treatment_category <- ifelse(grepl("DW", shannon_div_df$treatment), "Duckweed", "Water")
shannon_div_df$treatment_category <- factor(shannon_div_df$treatment_category, levels = c("Water", "Duckweed"))
# Mapping pond names to full names
pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson \nFarm Pond 1",
    "TF-2" = "Thompson \nFarm Pond 2",
    "UM-1" = "Upper \nMill Pond")

# Create a new pond_treatment column to control x-axis order
shannon_div_df <- shannon_div_df |> 
    mutate(pond_treatment = paste0(pond, "_", ifelse(treatment_category == "Water", "W", "D"))) |> 
    mutate(pond_treatment = factor(pond_treatment, levels = c(
        "MP-1_W", "MP-1_D",
        "UM-1_W", "UM-1_D",
        "ODR-2_W", "ODR-2_D",
        "ODR-3_W", "ODR-3_D",
        "TF-1_W", "TF-1_D",
        "TF-2_W", "TF-2_D"))) |> 
    mutate(pond_label = pond_name_mapping[pond])

# Calculate summary stats
div_stats <- shannon_div_df %>%
    group_by(pond, pond_label, treatment_category, pond_treatment) %>%
    summarise(mean_div = mean(ShannonDiversity),
              sd_div = sd(ShannonDiversity),
              .groups = "drop")

# Set dodge width
dodge_width <- 0.5

# Updated Plot using jitterdodge settings
div_plot <- ggplot(div_stats, aes(x = pond_label, y = mean_div, color = treatment_category)) +
    geom_jitter(data = shannon_div_df, 
                aes(x = pond_label, 
                    y = ShannonDiversity, 
                    color = treatment_category),
                size = 1.8, shape=1,
                position = position_jitterdodge(jitter.width = 0.2, 
                                                dodge.width = dodge_width),
                inherit.aes = FALSE) +
    geom_point(size = 2, shape = 19, 
               position = position_dodge(width = dodge_width)) +
    geom_errorbar(aes(ymin = mean_div - sd_div, ymax = mean_div + sd_div), 
                  width = 0.2, position = position_dodge(width = dodge_width)) +
    scale_color_manual(values = c("Water" = "#88CCEE", 
                                  "Duckweed" = "#44AA99"), 
                       name = "Sample Type") +
    labs(x = "Green Manure Source", 
         y = "Shannon Diversity Index") +
    scale_x_discrete(labels = unique(shannon_div_df$pond_label)) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))

div_plot

ggsave("Expt2_diversity_plot.jpg", div_plot, width = 6, height = 4)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt2_diversity_plot.jpg", div_plot, width = 6, height = 4)


#MCMCglmm and ANOVA Analysis Section ----------------------------------------

#MCMCglmm 
mcmc_data <- shannon_div_df %>%
    select(ShannonDiversity, pond, treatment_category) %>%
    mutate(pond = as.factor(pond),
           treatment_category = as.factor(treatment_category))

mod1 <- MCMCglmm(ShannonDiversity ~ treatment_category, data=mcmc_data,verbose=F, nitt = 101000, thin = 10, burnin = 1000) 
summary(mod1)

#are the ponds different between duckweed samples?

#filter df to duckweed only
mcmc_DW_data <- shannon_div_df %>%
    filter(treatment_category == "Duckweed") %>%
    select(ShannonDiversity, pond, treatment_category) %>%
    mutate(pond = as.factor(pond),
           treatment_category = as.factor(treatment_category))

mod2 <- MCMCglmm(ShannonDiversity ~ -1 + pond, data=mcmc_DW_data,verbose=F, nitt = 101000, thin = 10, burnin = 1000)
summary(mod2)

post_summ <- summary(mod2)$solutions
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

#filter df to water only
mcmc_W_data <- shannon_div_df %>%
    filter(treatment_category == "Water") %>%
    select(ShannonDiversity, pond, treatment_category) %>%
    mutate(pond = as.factor(pond),
           treatment_category = as.factor(treatment_category))

mod3 <- MCMCglmm(ShannonDiversity ~ -1 + pond, data=mcmc_W_data,verbose=F, nitt = 101000, thin = 10, burnin = 1000)
summary(mod3)

post_summ <- summary(mod3)$solutions
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

#average shannon score for duckweed treatment
mean(shannon_div_df$ShannonDiversity[shannon_div_df$treatment_category == "Duckweed"])
#average shannon score for water treatment
mean(shannon_div_df$ShannonDiversity[shannon_div_df$treatment_category == "Water"])

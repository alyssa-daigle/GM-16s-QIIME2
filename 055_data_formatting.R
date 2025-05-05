library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(vegan)
library(tibble)
library(MCMCglmm)

# Data Preparation Section --------------------------------------------------

#set working directory
setwd("/Users/alyssadaigle/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s")

# Step 1: Prepare OTU matrix
otu_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/feature_table_no_contam.tsv", header = TRUE, sep = "\t", check.names = FALSE)
colnames(otu_mat)[colnames(otu_mat) == "#OTU ID"] <- "otu"
otu_mat <- otu_mat |>  
    tibble::column_to_rownames("otu") |> 
    as.matrix()

# Write out processed OTU matrix
write.table(otu_mat, "processed_otu_matrix.tsv", quote = FALSE)
write.table(otu_mat, "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_otu_matrix.tsv", quote = FALSE)


# Step 2: Prepare taxonomy matrix
tax_mat <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/taxonomy_no_contam.tsv", header = TRUE, sep = "\t")

colnames(tax_mat)[colnames(tax_mat) == "Feature.ID"] <- "otu"
tax_mat <- tax_mat |> 
    tibble::column_to_rownames("otu")

# Split and clean up taxonomy
tax_mat <- tax_mat |> 
    separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") |> 
    mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus), ~sub("^.*?__(.*)", "\\1", .)))

tax_mat <- as.matrix(tax_mat)

# Write out processed taxonomy matrix
write.table(tax_mat, file = "processed_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA)
write.table(tax_mat, file = "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_taxonomy_matrix.tsv", sep = "\t", quote = FALSE, col.names = NA) 


# Step 3: Prepare sample data
samples_df <- read.table("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment1_16s/metadata.tsv", header = TRUE)
colnames(samples_df)[colnames(samples_df) == "sampleid"] <- "sample"
samples_df$sample <- rownames(samples_df)

# Write out processed sample metadata
write.table(samples_df, "processed_sample_metadata.tsv", quote = FALSE)
write.table(samples_df, "~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s/processed_sample_metadata.tsv", quote = FALSE)



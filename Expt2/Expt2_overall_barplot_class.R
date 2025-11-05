# Load libraries -----------------------------------------------------------

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(tibble)
library(pals)
library(Polychrome)

# Set working directory ----------------------------------------------------

setwd("/Users/alyssadaigle/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/Experiment2_16s")

# Load data ----------------------------------------------------------------

otu_mat <- read.table("processed_otu_matrix.tsv", header = TRUE, check.names = FALSE)
otu_mat <- as.matrix(otu_mat)

tax_mat <- read.table("processed_taxonomy_matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
tax_mat <- as.matrix(tax_mat)

samples_df <- read.table("processed_sample_metadata.tsv", header = TRUE)

# Filter OTU matrix --------------------------------------------------------

col_index <- which(colnames(otu_mat) == "MP-1-1")
otu_mat_subset <- otu_mat[, col_index:ncol(otu_mat)]
otu_mat_subset <- otu_mat_subset[, !grepl("Expt1|soil|blank", colnames(otu_mat_subset))]
otu_mat_subset_filtered <- otu_mat_subset[rowSums(otu_mat_subset) > 0, ]

# Filter taxonomy ----------------------------------------------------------

tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(otu_mat_subset_filtered), ]

# Initial filtering: remove only Mitochondria and fully unclassified families
tax_mat_filtered <- tax_mat[!(
    is.na(tax_mat[, "Family"]) |
        tax_mat[, "Class"] == "Mitochondria"
), ]

# Bin "Chloroplast" and "Incertae Sedis" into "Other"
tax_mat_filtered[, "Class"] <- ifelse(
    tax_mat_filtered[, "Class"] == "Chloroplast" | 
        grepl("Incertae Sedis", tax_mat_filtered[, "Class"], ignore.case = TRUE),
    "Other",
    tax_mat_filtered[, "Class"]
)


rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))
otu_mat_final <- otu_mat_subset_filtered[rownames(otu_mat_subset_filtered) %in% rownames(tax_mat_filtered), ]

# Filter samples -----------------------------------------------------------

row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[!grepl("Expt1|soil", rownames(samples_df_filtered)), , drop = FALSE]
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

# Reshape and merge --------------------------------------------------------

otu_mat_final_df <- as.data.frame(otu_mat_final) |> 
    rownames_to_column(var = "FeatureID")

otu_mat_long <- otu_mat_final_df |> 
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |> 
    left_join(samples_df_filtered, by = c("Sample" = "SampleID"))

tax_mat_filtered_df <- as.data.frame(tax_mat_filtered) |> 
    rownames_to_column(var = "FeatureID")

otu_mat_long <- otu_mat_long |> 
    left_join(tax_mat_filtered_df, by = "FeatureID")

# Filter out missing/unclassified ------------------------------------------

otu_mat_long_filtered <- otu_mat_long |> 
    filter(!is.na(Class) & Class != "Unclassified")

# Fill in missing combinations ---------------------------------------------

otu_treatment_counts_full <- otu_mat_long_filtered |> 
    complete(treatment, Class, fill = list(Count = 0))

otu_treatment_counts_full$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", otu_treatment_counts_full$treatment)
otu_treatment_counts_full$treatment_category <- ifelse(grepl("DW", otu_treatment_counts_full$treatment), "Duckweed", "Water")

otu_treatment_counts_full$treatment_category <- factor(otu_treatment_counts_full$treatment_category, levels = c("Water", "Duckweed"))

# Get top 15 classes -------------------------------------------------------

top15_classes <- otu_treatment_counts_full |> 
    group_by(Class) |> 
    summarise(TotalCount = sum(Count), .groups = "drop") |> 
    arrange(desc(TotalCount)) |> 
    slice_head(n = 16)

otu_treatment_counts_grouped <- otu_treatment_counts_full |> 
    mutate(Class = ifelse(Class %in% top15_classes$Class, Class, "Other")) |> 
    group_by(Class, pond, treatment_category) |> 
    summarise(Count = sum(Count), .groups = "drop")

# Relative abundance -------------------------------------------------------

otu_relative_abundance <- otu_treatment_counts_grouped |> 
    group_by(pond, treatment_category) |> 
    mutate(RelativeAbundance = Count / sum(Count)) |> 
    ungroup()

class_order <- otu_relative_abundance |> 
    group_by(Class) |> 
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |> 
    arrange(desc(TotalAbundance)) |> 
    pull(Class)


class_order <- as.character(class_order)
class_order <- class_order[class_order != "Other"]
class_order <- c("Other", rev(class_order))
otu_relative_abundance$Class <- factor(otu_relative_abundance$Class, levels = class_order, ordered = TRUE)

# Label mapping ------------------------------------------------------------

pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson Farm\nPond 1",
    "TF-2" = "Thompson Farm\nPond 2",
    "UM-1" = "Upper \nMill Pond")

otu_relative_abundance$pond_full_name <- recode(otu_relative_abundance$pond, !!!pond_name_mapping)

pond_levels <- c(
    "Mill Pond",
    "Upper \nMill Pond",
    "Dairy Farm\nPond 1",
    "Dairy Farm\nPond 2",
    "Thompson Farm\nPond 1",
    "Thompson Farm\nPond 2"
)

otu_relative_abundance$pond_full_name <- factor(otu_relative_abundance$pond_full_name, levels = pond_levels)

# Plotting -----------------------------------------------------------------

# Your palette
gravityFalls_colors <- c(
    "#474747FF",  # gray
    "#8B4513FF",  # brown
    "#D2B48CFF",  # tan (brownish)
    "#000000FF",  # black
    "#417BA1FF",  # blue
    "hotpink",  # pink
    "#FFFF2EFF",  # yellow
    "#345634FF",  # dark green
    "#8B0000FF",  # dark red
    "#E2725B",  # orange
    "#93C0D5FF",  # light blue
    "#9248A7FF",  # purple
    "#1C8859FF", # green
    "pink2", # pink
    "#8FBC8FFF" # light green
)


# Ensure 'Other' is first in family_order
class_order <- c("Other", setdiff(class_order, "Other"))

# Create final color vector: gray for "Other", rest from gravityFalls
custom_colors <- c("#999999", gravityFalls_colors[1:(length(class_order) - 1)])
names(custom_colors) <- class_order

p <- ggplot(otu_relative_abundance, aes(x = treatment_category, y = RelativeAbundance, fill = Class)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_cowplot() +
    labs(x = "Sample Type", y = "Relative Abundance", fill = "Bacterial Class") +
    facet_wrap(~pond_full_name, scales = "free_x", ncol = 6) + 
    scale_fill_manual(values = custom_colors, 
                      guide = guide_legend(ncol = 1, 
                                           keyheight = 0.75, 
                                           keywidth = 0.75,
                                           reverse = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          legend.position = "right",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"), 
          plot.margin = margin(10, 10, 10, 10)) 


print(p)

ggsave("Expt2_barplot_class.jpg", p, width = 14, height = 6)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt2_barplot_class.jpg", p, width = 14, height = 6)

# Compute average abundance by pond, treatment, and class, excluding "Other"
top4_by_pond_treatment <- otu_relative_abundance |>
    filter(Class != "Other") |>
    group_by(pond_full_name, treatment_category, Class) |>
    summarise(AvgAbundance = mean(RelativeAbundance), .groups = "drop") |>
    arrange(pond_full_name, treatment_category, desc(AvgAbundance)) |>
    group_by(pond_full_name, treatment_category) |>
    slice_head(n = 4) |>
    ungroup()

top4_by_pond_treatment |> 
    arrange(treatment_category, Class) |> 
    distinct(treatment_category, Class)

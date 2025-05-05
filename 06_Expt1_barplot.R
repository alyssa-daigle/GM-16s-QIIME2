library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(tibble)

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

# Filter metadata for Expt 1 specific analysis
row_index <- which(rownames(samples_df) == "24-d6")
samples_df_filtered <- samples_df[1:row_index, , drop = FALSE]

# Bar Plot Section -----------------------------------------------------------

# Step 1: Ensure missing families in treatments are treated as having zero counts

otu_mat_final_df <- as.data.frame(otu_mat_final) |> 
    rownames_to_column(var = "FeatureID")

# Step 2: Reshape OTU data into long format and merge with metadata
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

otu_mat_long <- otu_mat_final_df |> 
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |> 
    left_join(samples_df_filtered, by = c("Sample" = "SampleID"))

# Step 3: Add taxonomy information
tax_mat_filtered_df <- as.data.frame(tax_mat_filtered) |> 
    rownames_to_column(var = "FeatureID")

otu_mat_long <- otu_mat_long |> 
    left_join(tax_mat_filtered_df, by = "FeatureID")

# Step 4: Filter out unclassified taxa before calculating relative abundance

otu_mat_long_filtered <- otu_mat_long |> 
    filter(!is.na(Family) & Family != "Unclassified" & Family != "Mitochondria",
           !is.na(Class) & Class != "Unclassified")

# Step 5: Fill missing family counts with zeroes for each treatment
otu_treatment_counts_full <- otu_mat_long_filtered |> 
    complete(treatment, Family, fill = list(TotalCount = 0))  # Fill with zero for missing families

top15_families <- otu_treatment_counts_full |> 
    group_by(Family) |> 
    summarise(TotalCount = sum(Count), .groups = 'drop') |> 
    arrange(desc(TotalCount)) |> 
    slice_head(n = 15) |>
    ungroup()

full_top15_families <- otu_treatment_counts_full |> 
    filter(Family %in% top15_families$Family) 

#classify other families into "Other" group with total count
other_families <- otu_treatment_counts_full |> 
    mutate(Family = ifelse(Family %in% top15_families$Family, Family, "Other")) |> 
    filter(Family == "Other") |> 
    group_by(treatment, Family) |>
    summarise(Count = sum(Count), .groups = 'drop') 

#append "Other" to the top 15 DF, with Other coming as the last "Family" per treatment
full_top15_families <- full_top15_families |> 
    bind_rows(other_families) |> 
    mutate(Family = factor(Family, levels = c(top15_families$Family, "Other")))

#relative abundance per treatment
full_top15_families <- full_top15_families |> 
    group_by(treatment) |> 
    mutate(RelativeAbundance = Count / sum(Count)) |> 
    ungroup() |> 
    select(-FeatureID, -Kingdom, -Phylum, -Class, -Order, -Genus)

# Step 11: extract treatment components
full_top15_families <- full_top15_families |> 
    separate(treatment, into = c("geno", "cyano", "micro"), 
             sep = "-", remove = FALSE)

family_order <- full_top15_families |>
    group_by(Family) |>
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |>
    arrange(desc(TotalAbundance)) |>
    pull(Family)

# Convert to character to prevent NA coercion
family_order <- as.character(family_order)
family_order <- family_order[family_order != "Other"]
family_order <- c(rev(family_order))
family_order <- c("Other", family_order[family_order != "Other"])

print(family_order)
 
# Apply ordered factor to Family column
full_top15_families$Family <- factor(full_top15_families$Family, levels = family_order, ordered = TRUE)

# Assign gray to "Other", then generate the rest
glasbey_colors <- c("#999999", Polychrome::createPalette(length(family_order) - 1, seedcolors = c("#0000FF", "#00FF00", "#FF0000")))

# Name the colors according to the family order
names(glasbey_colors) <- family_order

full_top15_families$micro <- factor(full_top15_families$micro, levels = c("N", setdiff(unique(full_top15_families$micro), "N")))

micro_labels <- c("N" = "Uninoculated",
                  "H" = "Home microbiome",
                  "KF" = "Kingman Farm microbiome",
                  "ODR" = "Dairy Farm microbiome")

# Step 1: Define and preprocess the variables
# Create x_axis as a factor with fixed levels to preserve the exact order
full_top15_families$group <- sub("-[YN]-.*", "", full_top15_families$treatment)  # DR, LR, etc.
full_top15_families$inoc <- sub("^.*?-([YN])-.*", "\\1", full_top15_families$treatment)  # Y or N
full_top15_families$x_axis <- paste(full_top15_families$group, full_top15_families$inoc, sep = "_")

# Define levels for x_axis and labels for the x-axis
x_levels <- c("DR_Y", "DR_N", "LR_Y", "LR_N", "M_Y", "M_N", "TF_Y", "TF_N", "UM_Y", "UM_N", "W_Y", "W_N")
inoc_labels <- c("Y", "N", "Y", "N", "Y", "N", "Y", "N", "Y", "N", "Y", "N")

# Ensure the x_axis is a factor with the predefined levels
full_top15_families$x_axis <- factor(full_top15_families$x_axis, levels = x_levels)

# Set the labels for the x-axis, keeping Y/N only (no group labels like DR, LR, etc.)
x_labels <- setNames(inoc_labels, x_levels)

# Step 2: Create the plot
p <- ggplot(full_top15_families, aes(x = x_axis, y = RelativeAbundance, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_cowplot() +
    labs(x = "Treatment", y = "Relative Abundance", fill = "Family") +
    facet_wrap(~micro, scales = "free_x", ncol = 4, labeller = labeller(micro = micro_labels)) +
    scale_x_discrete(labels = x_labels) +
    theme(
        axis.text.x = element_text(size = 8, lineheight = 0.9),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "right",
        legend.text = element_text(size = 10)
    ) +
    scale_fill_manual(values = glasbey_colors, guide = guide_legend(reverse = TRUE, ncol = 1))

# Separate legend and plot
legend <- get_legend(p)
p_no_legend <- p + theme(legend.position = "none")
combined_plot <- plot_grid(
    p_no_legend, 
    NULL,  # adds a blank space
    legend, 
    ncol = 3, 
    rel_widths = c(0.82, 0.03, 0.15)
)

ggsave("Expt1_barplot.jpg", combined_plot, width = 14, height = 6)

ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_barplot.jpg", combined_plot, width = 14, height = 6)


library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(tibble)
library(pals)
library(Polychrome)

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


# Filter sample data for Expt 2 specific analysis
row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[!grepl("Expt1|soil", rownames(samples_df_filtered)), , drop = FALSE]

otu_mat_subset <- otu_mat_subset[, !grepl("Expt1|soil|blank", colnames(otu_mat_subset))]


# Bar Plot Section ---------------------------------------------------------

# Step 1: Ensure missing families in treatments are treated as having zero counts

otu_mat_final_df <- as.data.frame(otu_mat_final)|> 
    rownames_to_column(var = "FeatureID")

# Step 2: Reshape OTU data into long format and merge with metadata
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

otu_mat_long <- otu_mat_final_df|> 
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count")|> 
    left_join(samples_df_filtered, by = c("Sample" = "SampleID"))

# Step 3: Add taxonomy information
tax_mat_filtered_df <- as.data.frame(tax_mat_filtered)|> 
    rownames_to_column(var = "FeatureID")

otu_mat_long <- otu_mat_long|> 
    left_join(tax_mat_filtered_df, by = "FeatureID")


# Step 1: Remove unclassified and unwanted taxa
otu_mat_long_filtered <- otu_mat_long|> 
    filter(!is.na(Family) & Family != "Unclassified" & Family != "Mitochondria",
           !is.na(Class) & Class != "Unclassified")

# Step 2: Fill missing combinations
otu_treatment_counts_full <- otu_mat_long_filtered|> 
    complete(treatment, Family, fill = list(Count = 0))

# Extract pond and treatment category
otu_treatment_counts_full$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", otu_treatment_counts_full$treatment)
otu_treatment_counts_full$treatment_category <- ifelse(grepl("DW", otu_treatment_counts_full$treatment), "Duckweed", "Water")
otu_treatment_counts_full$treatment_category <- factor(otu_treatment_counts_full$treatment_category, levels = c("Water", "Duckweed"))


otu_cyano <- otu_treatment_counts_full %>%
    filter(Class == "Cyanophyceae")

# Step 3: Calculate top 15 families across all samples
top15_families <- otu_cyano|> 
    group_by(Family)|> 
    summarise(TotalCount = sum(Count), .groups = "drop")|> 
    arrange(desc(TotalCount))|> 
    slice_head(n = 15) #set to 16 to account for Insertae Sedis being moved to Other

# Step 4: Group other families as "Other"
otu_treatment_counts_grouped <- otu_cyano|> 
    mutate(Family = ifelse(Family %in% top15_families$Family & Family != "Incertae Sedis", Family, "Other")) |> 
    group_by(Family, pond, treatment_category)|> 
    summarise(Count = sum(Count), .groups = "drop")

# Step 5: Calculate relative abundance within each pond-treatment group
otu_relative_abundance <- otu_treatment_counts_grouped|> 
    group_by(pond, treatment_category)|> 
    mutate(RelativeAbundance = Count / sum(Count))|> 
    ungroup()

# Step 6: Order Family factor (Other first, then top 15 in reverse)
family_order <- otu_relative_abundance|> 
    group_by(Family)|> 
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop")|> 
    arrange(desc(TotalAbundance))|> 
    pull(Family)

family_order <- as.character(family_order)
family_order <- family_order[family_order != "Other"]
family_order <- c("Other", rev(family_order))
otu_relative_abundance$Family <- factor(otu_relative_abundance$Family, levels = family_order, ordered = TRUE)

# Mapping pond names to full names
pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson \nFarm Pond 1",
    "TF-2" = "Thompson \nFarm Pond 2",
    "UM-1" = "Upper \nMill Pond")

#create column pond_full_name
otu_relative_abundance$pond_full_name <- recode(otu_relative_abundance$pond, !!!pond_name_mapping)

pond_levels <- c(
    "MP-1" = "Mill Pond",
    "UM-1" = "Upper \nMill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson \nFarm Pond 1",
    "TF-2" = "Thompson \nFarm Pond 2")

otu_relative_abundance$pond_full_name <- factor(otu_relative_abundance$pond_full_name, levels = pond_levels)

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
class_order <- c("Other", setdiff(family_order, "Other"))

# Create final color vector: gray for "Other", rest from gravityFalls
custom_colors <- c("#999999", gravityFalls_colors[1:(length(family_order) - 1)])
names(custom_colors) <- family_order


# Step 9: Make the plot
p <- ggplot(otu_relative_abundance, aes(x = treatment_category, y = RelativeAbundance, fill = Family)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_cowplot() +
    labs(x = "Sample Type", y = "Relative Abundance", fill = "Cyanobacteria Family") +
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

# ggsave("Expt2_cyano_barplot.jpg", p, width = 14, height = 6)
# ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt2_cyano_barplot.jpg", p, width = 14, height = 6)


#which families were most abundant in water samples
otu_relative_abundance_water <- otu_relative_abundance|> 
    filter(treatment_category == "Water", Family != "Other")|> 
    group_by(Family)|> 
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop")|> 
    arrange(desc(TotalAbundance))

otu_relative_abundance_DW <- otu_relative_abundance|> 
    filter(treatment_category == "Duckweed",Family != "Other")|> 
    group_by(Family)|> 
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop")|> 
    arrange(desc(TotalAbundance))

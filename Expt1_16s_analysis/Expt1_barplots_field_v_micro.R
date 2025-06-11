library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(tibble)
library(scales)
library(stringr)
library(forcats)
library(ggtext)
library(MCMCglmm)


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
        grepl("W-Expt1", treatment)   ~ "Home",
        grepl("Dr-Expt1", treatment)  ~ "Home",
        grepl("LR-Expt1", treatment)  ~ "Home",
        grepl("UM-Expt1", treatment)  ~ "Home",
        grepl("MP-Expt1", treatment)  ~ "Home",
        grepl("TF-Expt1", treatment)  ~ "Home",
        grepl("KF-Expt1", treatment)  ~ "Kingman Farm",
        grepl("ODR-Expt1", treatment) ~ "Dairy Farm",
        TRUE ~ inoculum  # keep existing values for all other treatments
    ))

#rename treatment column
samples_df_filtered <- samples_df_filtered %>%
    mutate(treatment = case_when(
        grepl("Dr-Expt1", treatment) ~ gsub("Dr-Expt1", "A_DR-Expt1", treatment),
        grepl("LR-Expt1", treatment) ~ gsub("LR-Expt1", "A_LR-Expt1", treatment),
        grepl("MP-Expt1", treatment) ~ gsub("MP-Expt1", "A_MP-Expt1", treatment),
        grepl("TF-Expt1", treatment) ~ gsub("TF-Expt1", "A_TF-Expt1", treatment),
        grepl("UM-Expt1", treatment) ~ gsub("UM-Expt1", "A_UM-Expt1", treatment),
        grepl("W-Expt1", treatment) ~ gsub("W-Expt1", "A_W-Expt1", treatment), 
        grepl("KF-Expt1", treatment) ~ gsub("KF-Expt1", "A_KF-Expt1", treatment),
        grepl("ODR-Expt1", treatment) ~ gsub("ODR-Expt1", "A_ODR-Expt1", treatment),
        TRUE ~ treatment))

# Bar Plot Section -----------------------------------------------------------

# Step 1: Ensure missing families in treatments are treated as having zero counts

otu_mat_final_df <- as.data.frame(otu_mat_final) |> 
    rownames_to_column(var = "FeatureID")

# Step 2: Reshape OTU data into long format and merge with metadata
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

otu_mat_long <- otu_mat_final_df |> 
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |> 
    left_join(samples_df_filtered, by = c("Sample" = "SampleID")) |> 
    select(-Sample)

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
    complete(treatment, Family, fill = list(Count = 0))  # Fill missing counts with zero

# Step 6: Identify top 15 families by total count
top15_families <- otu_treatment_counts_full |>
    group_by(Family) |>
    summarise(TotalCount = sum(Count), .groups = 'drop') |>
    arrange(desc(TotalCount)) |>
    slice_head(n = 16)

# Step 7 & 8: Merge top 15 and "Other", and bin 'incertae sedis' into 'Other'
full_top15_families <- otu_treatment_counts_full |>
    mutate(Family = ifelse(Family %in% top15_families$Family & !grepl("incertae sedis", Family, ignore.case = TRUE),
                           Family,
                           "Other"))

# Add inoculum info (assumes treatment uniquely identifies inoculum)
inoculum_lookup <- samples_df_filtered |> 
    select(treatment, inoculum) |> 
    distinct()

full_top15_families <- full_top15_families |>
    left_join(inoculum_lookup, by = "treatment") |>
    select(-inoculum.x, inoculum = inoculum.y)


# Step 9: Summarize counts per sample, treatment, family, and inoculum
full_top15_families <- full_top15_families |>
    group_by(sample, treatment, Family, inoculum) |>
    summarise(Count = sum(Count), .groups = "drop") |> 
    ungroup()

# Step 10: Recalculate relative abundance per treatment
full_top15_families <- full_top15_families |>
    group_by(treatment) |>
    mutate(RelativeAbundance = Count / sum(Count)) |>
    ungroup()

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

#comparing "green manure farm" treatments
kf_odr_plot_df <- full_top15_families %>%
    filter(inoculum %in% c("Kingman Farm", "Dairy Farm"))

kf_odr_plot_df <- kf_odr_plot_df %>%
    mutate(
        treatment_label = case_when(
            str_starts(treatment, "A_") ~ "Field \nInoculum",
            TRUE ~ str_split(treatment, "-", simplify = TRUE)[, 2]
        ),
        treatment_label = factor(treatment_label, levels = c("Field \nInoculum", "Y", "N")),
        treatment_label = fct_recode(treatment_label,
                                     "+" = "Y",
                                     "-" = "N"
        )
    )

# Fix your named color vector (use = instead of ==)
colors <- c(
    "Weeksellaceae" = "#474747FF",       # gray
    "Caulobacteraceae" = "#8B4513FF",    # brown
    "Sphingobacteriaceae" = "#D2B48CFF", # tan
    "Methylophilaceae" = "#000000FF",    # black
    "Pseudobdellovibrionaceae" = "#417BA1FF", # blue
    "Chromobacteriaceae" = "hotpink",    # pink
    "Sphingomonadaceae" = "#FFFF2EFF",   # yellow
    "Spirosomataceae" = "#345634FF",     # dark green
    "Moraxellaceae" = "#8B0000FF",       # dark red
    "Dunaliellaceae" = "#E2725B",        # orange
    "Rhizobiaceae" = "#93C0D5FF",         # light blue
    "Flavobacteriaceae" = "#9248A7FF",   # purple
    "Enterobacteriaceae" = "#1C8859FF",  # green
    "Pseudomonadaceae" = "pink2",        # pink variant
    "Comamonadaceae" = "#8FBC8FFF",
    "Other" = "#999999"
)

# Plot using the combined palette
kf_odr_plot <- ggplot(kf_odr_plot_df, aes(x = treatment_label, 
                                          y = Count, 
                                          fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~inoculum, scales = "free_x") +
    labs(x = "Treatment", y = "Relative Abundance", fill = "Bacterial Family") +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
          axis.text.y = element_text(size = 8),
          legend.position = "right",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"), 
          plot.margin = margin(10, 10, 10, 10)) +
    scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE, ncol = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

kf_odr_plot

kf_odr_legend <- get_legend(kf_odr_plot)
kf_odr_plot_no_legend <- kf_odr_plot + theme(legend.position = "none")
kf_odr_combined_plot <- plot_grid(
    kf_odr_plot_no_legend, 
    kf_odr_legend, 
    ncol = 2, 
    rel_widths = c(0.65, 0.35)
)

kf_odr_combined_plot

ggsave("Expt1_kf_odr_combined_barplot.jpg", plot = kf_odr_combined_plot, width = 10, height = 6, dpi = 500)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_kf_odr_combined_barplot.jpg", plot = kf_odr_combined_plot, width = 10, height = 6, dpi = 500)




#comparing home treatments
home_field_plot_df <- full_top15_families %>%
    filter(grepl("Home", inoculum) | 
               (inoculum == "Home" & grepl("Expt1", treatment)))

home_field_plot_df <- home_field_plot_df |> 
    mutate(geno = case_when(
        grepl("(^|[^A-Za-z])[Dd]R\\b", treatment) ~ "DR",
        grepl("LR", treatment) ~ "LR",
        grepl("UM", treatment) ~ "UM",
        grepl("M", treatment) ~ "MP",
        grepl("TF", treatment) ~ "TF",
        grepl("W", treatment) ~ "W",
        TRUE ~ NA_character_
    ))

home_field_plot_df <- home_field_plot_df %>%
    mutate(
        treatment_label = case_when(
            str_starts(treatment, "A_") ~ "Field \nInoculum",
            TRUE ~ str_split(treatment, "-", simplify = TRUE)[, 2]
        ),
        treatment_label = factor(treatment_label, levels = c("Field \nInoculum", "Y", "N")),
        treatment_label = fct_recode(treatment_label,
                                     "+" = "Y",
                                     "-" = "N"
        )
    )

geno_labels <- c("DR" = "Durham \nReservoir", "LR" = "LaRoche \nPond", "MP" = "Mill \nPond", "UM" = "Upper \nMill Pond", "TF" = "Thompson \nFarm",  "W" = "Woodman \nRoad")

home_field_plot <- ggplot(home_field_plot_df, aes(x = treatment_label, y = Count, fill = Family)) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~geno, scales = "free_x", ncol = 6, labeller = labeller(geno = geno_labels)) +
    theme_cowplot() +
    labs(x = "Treatment", y = "Relative Abundance", fill = "Bacterial Family") +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, lineheight = 0.9),
          axis.text.y = element_text(size = 8),
          legend.position = "right",
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11, face = "bold"),
          strip.text = element_text(size = 10, face = "bold"), 
          plot.margin = margin(10, 10, 10, 10)) +
    scale_fill_manual(values = colors, guide = guide_legend(reverse = TRUE, ncol = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
home_field_plot

home_field_legend <- get_legend(home_field_plot)
home_field_no_legend <- home_field_plot + theme(legend.position = "none")
home_field_combined_plot <- plot_grid(
    home_field_no_legend, 
    # NULL,  # adds a blank space
    home_field_legend, 
    ncol = 3, 
    rel_widths = c(0.82, 0.03, 0.15)
)

home_field_combined_plot

ggsave("Expt1_home_field_combined_barplot.jpg", plot = home_field_combined_plot, width = 14, height = 6, dpi = 500)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_home_field_combined_barplot.jpg", plot = home_field_combined_plot, width = 14, height = 6, dpi = 500)

##Diversity analysis with MCMC -----------------------------------------

#Green manure farms first

#filter otu_mat_final to include only Kingman Farm and Dairy Farm inoculum

GM_shan_div <- full_top15_families %>%
    filter(inoculum %in% c("Kingman Farm", "Dairy Farm")) %>%
    group_by(treatment, inoculum) %>%
    summarise(Shannon = diversity(Count, index = "shannon"), .groups = 'drop')


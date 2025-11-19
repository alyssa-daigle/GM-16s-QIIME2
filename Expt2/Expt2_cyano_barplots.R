source("globals.R")
source("theme.R")

# Load .env file
load_dot_env()

# Get data path from environment variable
data_path <- Sys.getenv("data_path")
plots_path <- Sys.getenv("plots")

# Load data ----------------------------------------------------------------

# Read ASV matrix
asv_file <- file.path(data_path, "processed_asv_matrix.tsv")

asv_mat <- read.table(
    asv_file,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    check.names = FALSE
)
asv_mat <- as.matrix(asv_mat)

# Read taxonomy matrix
tax_file <- file.path(data_path, "processed_taxonomy_matrix.tsv")

tax_mat <- read.table(
    tax_file,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
tax_mat <- as.matrix(tax_mat)

# Read sample metadata
samples_file <- file.path(data_path, "processed_sample_metadata.tsv")

samples_df <- read.table(
    samples_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)
# Data Filtering Section ----------------------------------------------------

# Subset asv matrix after MP-1-1
col_index <- which(colnames(asv_mat) == "MP-1-1")
asv_mat_subset <- asv_mat[, col_index:ncol(asv_mat)]
asv_mat_subset <- asv_mat_subset[,
    !grepl("Expt1|soil|blank", colnames(asv_mat_subset))
]

# Remove rows where all values are 0
asv_mat_subset_filtered <- asv_mat_subset[rowSums(asv_mat_subset) > 0, ]

# Subset the taxonomy matrix
tax_mat_filtered <- tax_mat[
    rownames(tax_mat) %in% rownames(asv_mat_subset_filtered),
]

# Filter out specific taxa
tax_mat_filtered <- tax_mat_filtered[
    !((tax_mat_filtered[, "Kingdom"] == "Bacteria" &
        tax_mat_filtered[, "Phylum"] == "Cyanobacteriota" &
        tax_mat_filtered[, "Class"] == "Chloroplast" &
        tax_mat_filtered[, "Order"] == "Chloroplast" &
        is.na(tax_mat_filtered[, "Family"])) |
        is.na(tax_mat_filtered[, "Family"]) |
        tax_mat_filtered[, "Class"] == "Mitochondria"),
]

# Ensure unique taxa names
rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))

# Filter asv matrix based on taxonomy
asv_mat_final <- asv_mat_subset_filtered[
    rownames(asv_mat_subset_filtered) %in% rownames(tax_mat_filtered),
]


# Filter sample data for Expt 2 specific analysis
row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[
    !grepl("Expt1|soil", rownames(samples_df_filtered)),
    ,
    drop = FALSE
]

asv_mat_subset <- asv_mat_subset[,
    !grepl("Expt1|soil|blank", colnames(asv_mat_subset))
]

# Bar Plot Section ---------------------------------------------------------

# Ensure missing families in treatments are treated as having zero counts

asv_mat_final_df <- as.data.frame(asv_mat_final) |>
    rownames_to_column(var = "FeatureID")

# Reshape asv data into long format and merge with metadata
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

asv_mat_long <- asv_mat_final_df |>
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |>
    left_join(samples_df_filtered, by = c("Sample" = "SampleID"))

# Add taxonomy information
tax_mat_filtered_df <- as.data.frame(tax_mat_filtered) |>
    rownames_to_column(var = "FeatureID")

asv_mat_long <- asv_mat_long |>
    left_join(tax_mat_filtered_df, by = "FeatureID")


# Remove unclassified and unwanted taxa
asv_mat_long_filtered <- asv_mat_long |>
    filter(
        !is.na(Family) & Family != "Unclassified" & Family != "Mitochondria",
        !is.na(Class) & Class != "Unclassified"
    )

# Fill missing combinations
asv_treatment_counts_full <- asv_mat_long_filtered |>
    complete(treatment, Family, fill = list(Count = 0))

# Extract pond and treatment category
asv_treatment_counts_full$pond <- sub(
    "^([A-Za-z]+-\\d+)-.*$",
    "\\1",
    asv_treatment_counts_full$treatment
)
asv_treatment_counts_full$treatment_category <- ifelse(
    grepl("DW", asv_treatment_counts_full$treatment),
    "Duckweed",
    "Water"
)
asv_treatment_counts_full$treatment_category <- factor(
    asv_treatment_counts_full$treatment_category,
    levels = c("Water", "Duckweed")
)

asv_cyano <- asv_treatment_counts_full |>
    filter(Class == "Cyanophyceae")

# Calculate top 15 families across all samples
top15_families <- asv_cyano |>
    group_by(Family) |>
    summarise(TotalCount = sum(Count), .groups = "drop") |>
    arrange(desc(TotalCount)) |>
    slice_head(n = 15) #set to 16 to account for Insertae Sedis being moved to Other

# Group other families as "Other"
asv_treatment_counts_grouped <- asv_cyano |>
    mutate(
        Family = ifelse(
            Family %in% top15_families$Family & Family != "Incertae Sedis",
            Family,
            "Other"
        )
    ) |>
    group_by(Family, pond, treatment_category) |>
    summarise(Count = sum(Count), .groups = "drop")

# Calculate relative abundance within each pond-treatment group
asv_relative_abundance <- asv_treatment_counts_grouped |>
    group_by(pond, treatment_category) |>
    mutate(RelativeAbundance = Count / sum(Count)) |>
    ungroup()

# Order Family factor (Other first, then top 15 in reverse)
family_order <- asv_relative_abundance |>
    group_by(Family) |>
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |>
    arrange(desc(TotalAbundance)) |>
    pull(Family)

family_order <- as.character(family_order)
family_order <- family_order[family_order != "Other"]
family_order <- c("Other", rev(family_order))
asv_relative_abundance$Family <- factor(
    asv_relative_abundance$Family,
    levels = family_order,
    ordered = TRUE
)

#create column pond_full_name
asv_relative_abundance$pond_full_name <- recode(
    asv_relative_abundance$pond,
    !!!pond_name_mapping
)

asv_relative_abundance$pond_full_name <- factor(
    asv_relative_abundance$pond_full_name,
    levels = pond_levels
)

# Ensure 'Other' is first in family_order
class_order <- c("Other", setdiff(family_order, "Other"))

# Create final color vector: gray for "Other", rest from gravityFalls
custom_colors <- c("#999999", gravityFalls_colors[1:(length(family_order) - 1)])
names(custom_colors) <- family_order


# Make the plot
p <- ggplot(
    asv_relative_abundance,
    aes(x = treatment_category, y = RelativeAbundance, fill = Family)
) +
    geom_bar(stat = "identity", position = "stack") +
    theme_cowplot() +
    labs(
        x = "Sample Type",
        y = "Relative Abundance",
        fill = "Cyanobacteria Family"
    ) +
    facet_wrap(~pond_full_name, scales = "free_x", ncol = 6) +
    scale_fill_manual(
        values = custom_colors,
        guide = guide_legend(
            ncol = 1,
            keyheight = 0.75,
            keywidth = 0.75,
            reverse = TRUE
        )
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme(
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)
    )

print(p)

# Build full file path
barplot_plot_file <- file.path(plots_path, "Expt2_cyano_barplot.jpg")

# Save the plot
ggsave(
    filename = barplot_plot_file,
    plot = p,
    width = 14,
    height = 6,
    dpi = 500
)

#which families were most abundant in water samples
asv_relative_abundance_water <- asv_relative_abundance |>
    filter(treatment_category == "Water", Family != "Other") |>
    group_by(Family) |>
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |>
    arrange(desc(TotalAbundance))

asv_relative_abundance_DW <- asv_relative_abundance |>
    filter(treatment_category == "Duckweed", Family != "Other") |>
    group_by(Family) |>
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |>
    arrange(desc(TotalAbundance))

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

# Filter asv matrix --------------------------------------------------------

col_index <- which(colnames(asv_mat) == "MP-1-1")
asv_mat_subset <- asv_mat[, col_index:ncol(asv_mat)]
asv_mat_subset <- asv_mat_subset[,
    !grepl("Expt1|soil|blank", colnames(asv_mat_subset))
]
asv_mat_subset_filtered <- asv_mat_subset[rowSums(asv_mat_subset) > 0, ]

# Filter taxonomy ----------------------------------------------------------

tax_mat_filtered <- tax_mat[
    rownames(tax_mat) %in% rownames(asv_mat_subset_filtered),
]

# Initial filtering: remove only Mitochondria and fully unclassified families
tax_mat_filtered <- tax_mat[
    !(is.na(tax_mat[, "Family"]) |
        tax_mat[, "Class"] == "Mitochondria"),
]

# Bin "Chloroplast" and "Incertae Sedis" into "Other"
tax_mat_filtered[, "Class"] <- ifelse(
    tax_mat_filtered[, "Class"] == "Chloroplast" |
        grepl(
            "Incertae Sedis",
            tax_mat_filtered[, "Class"],
            ignore.case = TRUE
        ),
    "Other",
    tax_mat_filtered[, "Class"]
)

rownames(tax_mat_filtered) <- make.unique(rownames(tax_mat_filtered))
asv_mat_final <- asv_mat_subset_filtered[
    rownames(asv_mat_subset_filtered) %in% rownames(tax_mat_filtered),
]

# Filter samples -----------------------------------------------------------

row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[
    !grepl("Expt1|soil", rownames(samples_df_filtered)),
    ,
    drop = FALSE
]
samples_df_filtered$SampleID <- rownames(samples_df_filtered)

# Reshape and merge --------------------------------------------------------

asv_mat_final_df <- as.data.frame(asv_mat_final) |>
    rownames_to_column(var = "FeatureID")

asv_mat_long <- asv_mat_final_df |>
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |>
    left_join(samples_df_filtered, by = c("Sample" = "SampleID"))

tax_mat_filtered_df <- as.data.frame(tax_mat_filtered) |>
    rownames_to_column(var = "FeatureID")

asv_mat_long <- asv_mat_long |>
    left_join(tax_mat_filtered_df, by = "FeatureID")

# Filter out missing/unclassified ------------------------------------------

asv_mat_long_filtered <- asv_mat_long |>
    filter(!is.na(Class) & Class != "Unclassified")

# Fill in missing combinations ---------------------------------------------

asv_treatment_counts_full <- asv_mat_long_filtered |>
    complete(treatment, Class, fill = list(Count = 0))

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

# Get top 15 classes -------------------------------------------------------

top15_classes <- asv_treatment_counts_full |>
    group_by(Class) |>
    summarise(TotalCount = sum(Count), .groups = "drop") |>
    arrange(desc(TotalCount)) |>
    slice_head(n = 16)

asv_treatment_counts_grouped <- asv_treatment_counts_full |>
    mutate(Class = ifelse(Class %in% top15_classes$Class, Class, "Other")) |>
    group_by(Class, pond, treatment_category) |>
    summarise(Count = sum(Count), .groups = "drop")

# Relative abundance -------------------------------------------------------

asv_relative_abundance <- asv_treatment_counts_grouped |>
    group_by(pond, treatment_category) |>
    mutate(RelativeAbundance = Count / sum(Count)) |>
    ungroup()

class_order <- asv_relative_abundance |>
    group_by(Class) |>
    summarise(TotalAbundance = sum(RelativeAbundance), .groups = "drop") |>
    arrange(desc(TotalAbundance)) |>
    pull(Class)


class_order <- as.character(class_order)
class_order <- class_order[class_order != "Other"]
class_order <- c("Other", rev(class_order))
asv_relative_abundance$Class <- factor(
    asv_relative_abundance$Class,
    levels = class_order,
    ordered = TRUE
)

# Label mapping ------------------------------------------------------------

asv_relative_abundance$pond_full_name <- recode(
    asv_relative_abundance$pond,
    !!!pond_name_mapping
)

asv_relative_abundance$pond_full_name <- factor(
    asv_relative_abundance$pond_full_name,
    levels = pond_levels
)

# Plotting -----------------------------------------------------------------

# Ensure 'Other' is first in family_order
class_order <- c("Other", setdiff(class_order, "Other"))

# Create final color vector: gray for "Other", rest from gravityFalls
custom_colors <- c("#999999", gravityFalls_colors[1:(length(class_order) - 1)])
names(custom_colors) <- class_order

p <- ggplot(
    asv_relative_abundance,
    aes(x = treatment_category, y = RelativeAbundance, fill = Class)
) +
    geom_bar(stat = "identity", position = "stack") +
    theme_cowplot() +
    labs(
        x = "Sample Type",
        y = "Relative Abundance",
        fill = "Bacterial Class"
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
barplot_plot_file <- file.path(plots_path, "Expt2_barplot_class.jpg")

# Save the plot
ggsave(
    filename = barplot_plot_file,
    plot = p,
    width = 14,
    height = 6,
    dpi = 500
)

# Compute average abundance by pond, treatment, and class, excluding "Other"
top4_by_pond_treatment <- asv_relative_abundance |>
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

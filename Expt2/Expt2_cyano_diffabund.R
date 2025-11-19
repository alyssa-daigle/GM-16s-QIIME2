source("globals.R")
source("theme.R")

# Load .env file
load_dot_env()

# Get data path from environment variable
data_path <- Sys.getenv("data_path")
plots_path <- Sys.getenv("plots")

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

# Subset asv matrix after MP-1-1 (expt2 data only)
col_index <- which(colnames(asv_mat) == "MP-1-1")
asv_mat_subset <- asv_mat[, col_index:ncol(asv_mat)]
asv_mat_subset <- asv_mat_subset[,
    !grepl("Expt1|soil|blank", colnames(asv_mat_subset))
]
asv_mat_subset <- asv_mat_subset[rowSums(asv_mat_subset) > 0, ]

# Filter taxonomy to match asv
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(asv_mat_subset), ]

# Filter for Cyanophyceae
tax_mat_cyano <- tax_mat_filtered |>
    as.data.frame() |>
    filter(Class == "Cyanophyceae") |>
    drop_na(Family) |>
    filter(Family != c("Unclassified", "Unassigned"))

# Keep only rows in asv table that match Cyanophyceae (cyano-specific analysis)
asv_mat_cyano <- asv_mat_subset[
    rownames(asv_mat_subset) %in% rownames(tax_mat_cyano),
]

# Convert to long format and merge taxonomy early
asv_long_cyano <- as.data.frame(asv_mat_cyano) |>
    rownames_to_column("FeatureID") |>
    pivot_longer(-FeatureID, names_to = "Sample", values_to = "Count") |>
    left_join(
        tax_mat_cyano |> rownames_to_column("FeatureID"),
        by = "FeatureID"
    )

# Add metadata (filtered only once)
samples_df_filtered <- samples_df |>
    (\(df) df[!grepl("Expt1|soil|blank", rownames(df)), ])() |>
    (\(df) transform(df, SampleID = rownames(df)))()

# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
asv_mat_relative <- sweep(asv_mat_cyano, 2, colSums(asv_mat_cyano), `/`) * 100

# Filter taxa based on mean relative abundance
keep_taxa <- rownames(asv_mat_relative)[apply(asv_mat_relative, 1, function(x) {
    any(x >= 1)
})]
asv_mat_filtered_counts <- asv_mat_cyano[keep_taxa, ]

# Filter taxonomy again based on newly calculated rel abundance in the asv matrix
tax_mat_filtered_counts <- tax_mat_cyano[
    rownames(tax_mat_cyano) %in% rownames(asv_mat_filtered_counts),
]

# Phyloseq Object Creation --------------------------------------------------

# Create phyloseq object
ASV <- otu_table(asv_mat_filtered_counts, taxa_are_rows = TRUE)
tax_mat_filtered_counts <- as.matrix(tax_mat_filtered_counts)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(ASV, TAX, samples)

# Extract sample data
sample_data_df <- data.frame(phyloseq::sample_data(physeq))

# Create new column based on treatment values
sample_data_df$treatment_type <- factor(
    ifelse(grepl("DW$", sample_data_df$treatment), "duckweed", "water"),
    levels = c("water", "duckweed")
)

# Update the phyloseq object
phyloseq::sample_data(physeq) <- sample_data_df

physeq_family <- tax_glom(physeq, taxrank = "Family")

# Run the corncob test on the filtered phyloseq object
diff_test <- differentialTest(
    formula = ~treatment_type,
    ,
    phi.formula = ~treatment_type,
    formula_null = ~1,
    phi.formula_null = ~treatment_type,
    test = "Wald",
    boot = FALSE,
    data = physeq_family,
    fdr_cutoff = 0.05
)

# Extract the plot data
plot_data <- plot(diff_test, level = "Family", data_only = TRUE)

# Create a new column for enrichment direction
plot_data$effect_direction <- ifelse(plot_data$x > 0, "Enriched", "Depleted")

# Optional: reorder taxa for better plotting (by effect size within each facet)
plot_data <- plot_data |>
    group_by(variable) |>
    mutate(taxa = forcats::fct_reorder(taxa, x)) |>
    ungroup()

Expt2_cyano_diffabund_plot <- ggplot(
    plot_data,
    aes(x = x, y = taxa, xmin = xmin, xmax = xmax, color = effect_direction)
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(position = position_dodge(width = 0.0), size = 0.3) +
    scale_color_manual(values = c("Enriched" = "blue", "Depleted" = "red")) +
    theme_minimal(base_size = 12) +
    labs(
        x = "log(Odds Ratio) relative to Water Samples",
        y = "Cyanobacterial Family",
        color = "Direction of Effect"
    ) +
    theme_bw() +
    theme(
        axis.text.y = element_text(
            size = 8,
            margin = margin(l = 0.4, r = 0.1, unit = "cm")
        ),
        axis.text.x = element_text(angle = 0, hjust = 0.5)
    )

# Build full file path
diffabund_plot_file <- file.path(plots_path, "Expt2_cyano_diffabund_plot.jpg")

# Save the plot
ggsave(
    filename = diffabund_plot_file,
    plot = Expt2_cyano_diffabund_plot,
    width = 5,
    height = 4,
    dpi = 500
)

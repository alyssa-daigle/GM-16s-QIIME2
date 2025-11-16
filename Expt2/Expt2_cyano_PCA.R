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

# PCA Analysis ------------------------------------------------

asv_data <- otu_table(physeq)
asv_data_log <- log1p(asv_data)

# Perform PCA
pca_results <- rda(t(asv_data_log)) # Transpose the data for PCA
eig_vals <- eigenvals(pca_results) # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100 # Variance explained

# Extract PCA scores
pca_scores <- as.data.frame(scores(pca_results, display = "sites"))
pca_scores$SampleID <- rownames(pca_scores)
pca_scores <- merge(
    pca_scores,
    samples_df_filtered,
    by.x = "SampleID",
    by.y = "row.names",
    all.x = TRUE
)

# Extract pond category (first part before second dash)
pca_scores$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", pca_scores$treatment)

# Categorize treatment as 'Water' or 'Duckweed'
pca_scores$treatment_category <- ifelse(
    grepl("DW", pca_scores$treatment),
    "Duckweed",
    "Water"
)

# Recode pond names using pond_name_mapping
pca_scores$pond_full_name <- recode(pca_scores$pond, !!!pond_name_mapping)

# Apply the selected colors in ggplot
pca_plot <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = pond_full_name, shape = treatment_category)
) +
    geom_point(size = 3.5) +
    scale_color_manual(values = expt2_custom_colors) +
    labs(
        x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
        y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = ""),
        color = "Green Manure Source",
        shape = "Sample Type"
    ) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme_cowplot() +
    pca_theme()

print(pca_plot)

# Build full file path
pca_plot_file <- file.path(plots_path, "Expt2_cyano_pca_plot.jpg")

# Save the plot
ggsave(
    filename = pca_plot_file,
    plot = pca_plot,
    width = 6,
    height = 4,
    dpi = 500
)

# Duckweed-specific PCA plot --------------------------

# filter the physeq object to just duckweed samples and perform new PCA
physeq_dw <- subset_samples(physeq, grepl("DW", treatment))
# Perform PCA on duckweed samples
dw_asv_data <- otu_table(physeq_dw)
dw_asv_data_log <- log1p(dw_asv_data)

# Step 2: Perform PCA
pca_results <- rda(t(dw_asv_data_log)) # Transpose the data for PCA
eig_vals <- eigenvals(pca_results) # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100 # Variance explained

# Step 3: Extract PCA scores
pca_scores <- as.data.frame(scores(pca_results, display = "sites"))
pca_scores$SampleID <- rownames(pca_scores)
pca_scores <- merge(
    pca_scores,
    samples_df_filtered,
    by.x = "SampleID",
    by.y = "row.names",
    all.x = TRUE
)

# Mapping pond names to full names
pond_name_mapping <- c(
    "MP-1" = "Mill Pond",
    "ODR-2" = "Dairy Farm\nPond 1",
    "ODR-3" = "Dairy Farm\nPond 2",
    "TF-1" = "Thompson Farm\nPond 1",
    "TF-2" = "Thompson Farm\nPond 2",
    "UM-1" = "Upper Mill Pond"
)

#separate treatment into pond
pca_scores$pond <- sub("^([A-Za-z]+-\\d+)-.*$", "\\1", pca_scores$treatment)

pca_scores$pond <- recode(pca_scores$pond, !!!pond_name_mapping)

# Apply the selected colors in ggplot
dw_pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = pond)) +
    geom_point(size = 3.5) +
    scale_color_manual(
        values = expt2_custom_colors,
        labels = pond_name_mapping
    ) + # Use manually selected colors
    labs(
        x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
        y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = ""),
        color = "Green Manure Source"
    ) +
    theme_cowplot() +
    pca_theme()

# -----------------------------------------------------------
# which taxa are driving the PCA
# also linear model to see that PC1 is significantly affected by treatmemnt

# Load species (asv) scores
asv_loadings <- scores(pca_results, display = "species") |>
    as.data.frame() |>
    tibble::rownames_to_column("asv")

# Extract taxonomy table
taxonomy_df <- tax_table(physeq) |>
    as.data.frame() |>
    tibble::rownames_to_column("asv") |>
    select(asv, Family)

# Function to get top N asvs for a given PC
get_top_asvs <- function(df, pc_col, n = 10) {
    df |>
        arrange(desc(abs(.data[[pc_col]]))) |>
        slice_head(n = n) |>
        left_join(taxonomy_df, by = "asv") |>
        select(asv, !!pc_col, Family)
}

# Get top 10 for PC1 and PC2
top10_PC1_df <- get_top_asvs(asv_loadings, "PC1")
top10_PC2_df <- get_top_asvs(asv_loadings, "PC2")

# models -----------------------------------

# Ensure treatment is a factor (if not already)
pca_scores$treatment_category <- as.factor(pca_scores$treatment_category)

PC1 <- MCMCglmm(
    PC1 ~ -1 + treatment_category:pond,
    data = pca_scores,
    nitt = 101000,
    burnin = 1000,
    thin = 10,
    verbose = FALSE
)
summary(PC1)

PC2 <- MCMCglmm(
    PC2 ~ treatment_category:pond,
    data = pca_scores,
    nitt = 13000,
    burnin = 3000,
    thin = 10,
    verbose = FALSE
)
summary(PC2)

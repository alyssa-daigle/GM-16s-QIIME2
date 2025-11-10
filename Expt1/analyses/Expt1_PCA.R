source("globals.R")

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

# Subset ASV matrix up to "24-d6"
col_index <- which(colnames(asv_mat) == "24-d6")
asv_mat_subset <- asv_mat[, 1:col_index]

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

# Filter ASV matrix based on taxonomy
asv_mat_final <- asv_mat_subset_filtered[
    rownames(asv_mat_subset_filtered) %in% rownames(tax_mat_filtered),
]

# Filter sample data for Expt 1 specific analysis
row_index <- which(rownames(samples_df) == "24-d6")
samples_df_filtered <- samples_df[1:row_index, , drop = FALSE]


# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
asv_mat_relative <- sweep(asv_mat_final, 2, colSums(asv_mat_final), `/`) * 100

# Filter taxa based on mean relative abundance (keep only those with >= 1%)
keep_taxa <- rownames(asv_mat_relative)[apply(asv_mat_relative, 1, function(x) {
    any(x >= 1)
})]
asv_mat_filtered_counts <- asv_mat_final[keep_taxa, ]

# Filter taxonomy again based on newly calculated rel abundance in the ASV matrix
tax_mat_filtered_counts <- tax_mat_filtered[
    rownames(tax_mat_filtered) %in% rownames(asv_mat_filtered_counts),
]

# Phyloseq Object Creation --------------------------------------------------

# Create phyloseq object
ASV <- otu_table(asv_mat_filtered_counts, taxa_are_rows = TRUE) #actually ASVs
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(ASV, TAX, samples)

# PCA Analysis Section -------------------------------------------------------
asv_data <- otu_table(physeq)
asv_data_log <- log1p(asv_data)

# Perform PCA using prcomp()
pca_results <- prcomp(t(asv_data_log), scale. = TRUE)
eig_vals <- eigenvals(pca_results) # Eigenvalues
varexpl <- (eig_vals / sum(eig_vals)) * 100 # Variance explained

# Extract PCA scores for the first three components (PC1, PC2, and PC3)
pca_scores <- as.data.frame(pca_results$x[, 1:3])

# merge with your sample metadata:
pca_scores$SampleID <- rownames(pca_scores)

pca_scores <- merge(
    pca_scores,
    samples_df_filtered,
    by.x = "SampleID",
    by.y = "row.names",
    all.x = TRUE
)

# separate treatment into components
pca_scores <- pca_scores |>
    separate(
        treatment,
        into = c("geno", "cyano", "micro"),
        sep = "-",
        remove = FALSE
    )

# assign colors to microbes
custom_colors <- c("#AA4499", "#DDCC77", "#88CCEE", "#117733")

# Visualize PCA
pca_plot <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = micro, shape = cyano)
) +
    geom_point(size = 2.5) +
    scale_color_manual(
        name = "Microbiome Source",
        values = custom_colors,
        breaks = c("H", "KF", "ODR", "N"),
        labels = c("Home", "Kingman Farm", "Dairy Farm", "Uninoculated")
    ) +
    scale_shape_manual(
        name = "*M. aeruginosa* Spike",
        values = c("Y" = 16, "N" = 17),
        labels = c("Y" = "Yes", "N" = "No")
    ) +
    labs(
        x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
        y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = "")
    ) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme_cowplot() +
    theme(
        legend.title = element_markdown(size = 10),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 5.7)
    )

pca_plot

# Build full file path
pca_plot_file <- file.path(plots_path, "Expt1_pca_plot.jpg")

# Save the plot
ggsave(
    filename = pca_plot_file,
    plot = pca_plot,
    width = 6,
    height = 4,
    dpi = 500
)

# which taxa are driving the PCA

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


# linear model to see that PC1 is significantly affected by treatment

pca_scores$treatment <- as.factor(pca_scores$treatment)

lm2_PC1 <- MCMCglmm(
    PC1 ~ -1 + cyano + micro + geno,
    data = pca_scores,
    nitt = 101000,
    burnin = 1000,
    thin = 10,
    verbose = FALSE
)
summary(lm2_PC1)

lm_PC2 <- MCMCglmm(
    PC2 ~ -1 + cyano + micro + geno,
    data = pca_scores,
    nitt = 101000,
    burnin = 1000,
    thin = 10,
    verbose = FALSE
)
summary(lm_PC2)

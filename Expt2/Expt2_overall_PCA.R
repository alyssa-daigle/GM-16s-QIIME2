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

# Subset ASV matrix after MP-1-1
col_index <- which(colnames(asv_mat) == "MP-1-1")
asv_mat_subset <- asv_mat[, col_index:ncol(asv_mat)]
asv_mat_subset <- asv_mat_subset[,
    !grepl("Expt1|soil|blank", colnames(asv_mat_subset))
]

#  Remove rows where all values are 0
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

# Filter sample data for Expt 2 specific analysis
row_index <- which(rownames(samples_df) == "ODR-2-1")
samples_df_filtered <- samples_df[row_index:nrow(samples_df), , drop = FALSE]
samples_df_filtered <- samples_df_filtered[
    !grepl("Expt1|soil|blank", rownames(samples_df_filtered)),
    ,
    drop = FALSE
]

# Relative Abundance Calculation ----------------------------------------

# Relative abundance for filtering
asv_mat_relative <- sweep(asv_mat_final, 2, colSums(asv_mat_final), `/`) * 100

keep_taxa <- rownames(asv_mat_relative)[apply(asv_mat_relative, 1, function(x) {
    any(x >= 1)
})]
asv_mat_filtered_counts <- asv_mat_final[keep_taxa, ]

# Filter taxonomy again based on newly calculated rel abundance in the asv matrix
tax_mat_filtered_counts <- tax_mat_filtered[
    rownames(tax_mat_filtered) %in% rownames(asv_mat_filtered_counts),
]

# Phyloseq Object Creation --------------------------------------------------
# Create phyloseq object
ASV <- otu_table(asv_mat_filtered_counts, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat_filtered_counts)
samples <- sample_data(samples_df_filtered)
physeq <- phyloseq(ASV, TAX, samples)

# PCA Analysis Section -------------------------------------------------------

# Extract asv table and transform data
asv_data <- otu_table(physeq)
asv_data_log <- log1p(asv_data) # log(1 + count) transformation

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
    geom_point(size = 3.2) +
    scale_color_manual(values = expt2_custom_colors) + # Use manually selected colors
    labs(
        x = paste("PC1 (", round(varexpl[1], 2), "%)", sep = ""),
        y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = ""),
        color = "Green Manure \nSource",
        shape = "Sample Type"
    ) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
    guides(color = guide_legend(order = 1), shape = guide_legend(order = 2)) +
    theme_cowplot() +
    pca_theme())

print(pca_plot)

# Build full file path
pca_plot_file <- file.path(plots_path, "Expt2_pca_plot.jpg")

# Save the plot
ggsave(
    filename = pca_plot_file,
    plot = pca_plot,
    width = 6,
    height = 4,
    dpi = 500
)

###PC1 by microbiome source

# Calculate mean and standard deviation for each group of pond and sample type
pca_stats <- pca_scores |>
    group_by(pond, treatment_category) |>
    summarise(mean_PC1 = mean(PC1), sd_PC1 = sd(PC1), .groups = 'drop')

# Create the plot
# Set custom factor order for pond
custom_pond_order <- c("MP-1", "UM-1", "ODR-2", "ODR-3", "TF-1", "TF-2")
pca_stats$pond <- factor(pca_stats$pond, levels = custom_pond_order)
pca_scores$pond <- factor(pca_scores$pond, levels = custom_pond_order)

# Define the dodge width to reuse consistently
dodge_width <- 0.5

#plot
pc1_pond_plot <- ggplot(
    pca_stats,
    aes(x = pond, y = mean_PC1, color = treatment_category)
) +
    geom_jitter(
        data = pca_scores,
        aes(x = pond, y = PC1, color = treatment_category),
        size = 1.8,
        shape = 1,
        position = position_jitterdodge(
            jitter.width = 0.2,
            dodge.width = dodge_width
        ),
        inherit.aes = FALSE
    ) +
    geom_point(
        size = 2,
        shape = 19,
        position = position_dodge(width = dodge_width)
    ) +
    geom_errorbar(
        aes(ymin = mean_PC1 - sd_PC1, ymax = mean_PC1 + sd_PC1),
        width = 0.2,
        position = position_dodge(width = dodge_width)
    ) +
    scale_color_manual(
        name = "Sample Type",
        values = c("Water" = "#88CCEE", "Duckweed" = "#44AA99")
    ) +
    scale_x_discrete(labels = pond_name_mapping) +
    labs(
        x = "Green Manure Source",
        y = paste("PC1 (", round(varexpl[1], 2), "%)", sep = "")
    ) +
    theme_cowplot() +
    theme(
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
    )

pc1_pond_plot

###PC2 by microbiome source

# Calculate mean and standard deviation for PC2
pca_stats_pc2 <- pca_scores |>
    group_by(pond, treatment_category) |>
    summarise(mean_PC2 = mean(PC2), sd_PC2 = sd(PC2), .groups = 'drop')

# Set custom pond order again for PC2 data
pca_stats_pc2$pond <- factor(pca_stats_pc2$pond, levels = custom_pond_order)
pca_scores$pond <- factor(pca_scores$pond, levels = custom_pond_order) # This is already done but reassigning is safe

# Plot PC2
pc2_pond_plot <- ggplot(
    pca_stats_pc2,
    aes(x = pond, y = mean_PC2, color = treatment_category)
) +
    geom_jitter(
        data = pca_scores,
        aes(x = pond, y = PC2, color = treatment_category),
        size = 1.8,
        shape = 1,
        position = position_jitterdodge(
            jitter.width = 0.2,
            dodge.width = dodge_width
        ),
        inherit.aes = FALSE
    ) +
    geom_point(
        size = 2,
        shape = 19,
        position = position_dodge(width = dodge_width)
    ) +
    geom_errorbar(
        aes(ymin = mean_PC2 - sd_PC2, ymax = mean_PC2 + sd_PC2),
        width = 0.2,
        position = position_dodge(width = dodge_width)
    ) +
    scale_color_manual(
        name = "Sample Type",
        values = c("Water" = "#88CCEE", "Duckweed" = "#44AA99")
    ) +
    scale_x_discrete(labels = pond_name_mapping) +
    labs(
        x = "Green Manure Source",
        y = paste("PC2 (", round(varexpl[2], 2), "%)", sep = "")
    ) +
    theme_cowplot() +
    theme(
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
    )

# Display the plot
pc2_pond_plot

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


# Ensure treatment is a factor (if not already)
pca_scores$treatment_category <- as.factor(pca_scores$treatment_category)

PC1 <- MCMCglmm(
    PC1 ~ -1 + treatment_category:pond,
    data = pca_scores,
    nitt = 13000,
    burnin = 3000,
    thin = 10,
    verbose = FALSE
)
summary(PC1)

PC2 <- MCMCglmm(
    PC2 ~ -1 + treatment_category:pond,
    data = pca_scores,
    nitt = 13000,
    burnin = 3000,
    thin = 10,
    verbose = FALSE
)
summary(PC1)

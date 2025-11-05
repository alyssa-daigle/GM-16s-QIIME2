source("globals.R")
load_dot_env()

# Get data path
data_path <- Sys.getenv("data_path")
plots_path <- Sys.getenv("plots")

# -------------------------
# Read ASV matrix
# -------------------------
asv_file <- file.path(data_path, "processed_asv_matrix.tsv")
asv_mat <- read.table(
  asv_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
asv_mat <- as.matrix(asv_mat)

# -------------------------
# Read taxonomy matrix
# -------------------------
tax_file <- file.path(data_path, "processed_taxonomy_matrix.tsv")
tax_mat <- read.table(
  tax_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)
tax_mat <- as.matrix(tax_mat)

# -------------------------
# Read sample metadata
# -------------------------
samples_file <- file.path(data_path, "processed_sample_metadata.tsv")
samples_df <- read.table(
  samples_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# -------------------------
# Filter samples (Expt1, up to "24-d6")
# -------------------------
col_index <- which(colnames(asv_mat) == "24-d6")
asv_mat_subset <- asv_mat[, 1:col_index]

# Remove ASVs with all zeros
asv_mat_subset <- asv_mat_subset[rowSums(asv_mat_subset) > 0, ]

# Filter taxonomy to keep only ASVs in filtered ASV matrix
tax_mat_filtered <- tax_mat[rownames(tax_mat) %in% rownames(asv_mat_subset), ]

# -------------------------
# Remove Uninoculated samples from metadata
# -------------------------
samples_df_filtered <- samples_df[1:col_index, ] |>
  separate(
    treatment,
    into = c("geno", "cyano", "micro"),
    sep = "-",
    remove = FALSE
  ) |>
  filter(micro != "N") # Remove Uninoculated

# Keep only the corresponding columns in ASV matrix
asv_mat_filtered_samples <- asv_mat_subset[, samples_df_filtered$sample]

# -------------------------
# CLR Transformation
# -------------------------
# Add pseudocount to avoid log(0)
asv_mat_clr <- t(apply(asv_mat_filtered_samples + 1, 2, function(x) clr(x)))

# -------------------------
# PCA
# -------------------------
pca_results <- prcomp(asv_mat_clr, scale. = FALSE)

# Variance explained
varexpl <- (pca_results$sdev^2 / sum(pca_results$sdev^2)) * 100

# PCA scores for PC1 and PC2
pca_scores <- as.data.frame(pca_results$x[, 1:3])
pca_scores$SampleID <- rownames(pca_scores)

# Merge with sample metadata
pca_scores <- pca_scores |>
  left_join(samples_df_filtered, by = c("SampleID" = "sample"))

# -------------------------
# Plot PCA
# -------------------------
custom_colors <- c("#AA4499", "#DDCC77", "#88CCEE", "#117733") # Microbiome sources

pca_plot_noUninoc <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = micro, shape = cyano)
) +
  geom_point(size = 3) +
  scale_color_manual(
    name = "Microbiome Source",
    values = custom_colors,
    breaks = c("H", "KF", "ODR"),
    labels = c("Home", "Kingman Farm", "Dairy Farm")
  ) +
  scale_shape_manual(
    name = "*M. aeruginosa* Spike",
    values = c("Y" = 16, "N" = 17),
    labels = c("Y" = "Yes", "N" = "No")
  ) +
  labs(
    x = paste0("PC1 (", round(varexpl[1], 2), "%)"),
    y = paste0("PC2 (", round(varexpl[2], 2), "%)")
  ) +
  theme_cowplot() +
  theme(
    legend.title = element_markdown(size = 10),
    legend.text = element_text(size = 8)
  )

pca_plot_noUninoc

# Build full file path
pca_plot_file <- file.path(plots_path, "Expt1_pca_plot_noUninoc.jpg")

# Save the plot
ggsave(
  filename = pca_plot_file,
  plot = pca_plot_noUninoc,
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

#  linear model to see that PC1 is significantly affected by treatment

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

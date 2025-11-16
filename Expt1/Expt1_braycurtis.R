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

samples_df_filtered <- samples_df_filtered |>
  separate(
    treatment,
    into = c("geno", "cyano", "micro"),
    sep = "-",
    remove = FALSE
  )

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

# Bray-Curtis Dissimilarity Analysis --------------------------------------------------
bray_dist <- phyloseq::distance(physeq, method = "bray")

ordu <- ordinate(physeq, method = "PCoA", distance = bray_dist)

plot_ordination(physeq, ordu, color = "micro", shape = "cyano") +
  geom_point(size = 1) +
  theme_bw()


# Convert Bray–Curtis distance to long format
bray_df <- as.data.frame(as.matrix(bray_dist))
bray_df$Sample1 <- rownames(bray_df)

bray_long <- bray_df |>
  pivot_longer(
    cols = -Sample1,
    names_to = "Sample2",
    values_to = "BrayCurtis"
  ) |>
  filter(Sample1 != Sample2)

# Join metadata (geno, cyano, micro)
meta <- samples_df_filtered
meta$SampleID <- rownames(meta)

bray_long2 <- bray_long |>
  left_join(meta, by = c("Sample1" = "SampleID")) |>
  rename(
    geno1 = geno,
    cyano1 = cyano,
    micro1 = micro
  ) |>
  left_join(meta, by = c("Sample2" = "SampleID")) |>
  rename(
    geno2 = geno,
    cyano2 = cyano,
    micro2 = micro
  )

# Define Within vs Between comparisons for each factor
bray_long2 <- bray_long2 |>
  mutate(
    micro_pair = ifelse(micro1 == micro2, "Within-micro", "Between-micro"),
    cyano_pair = ifelse(cyano1 == cyano2, "Within-cyano", "Between-cyano"),
    geno_pair = ifelse(geno1 == geno2, "Within-geno", "Between-geno")
  )

# Convert to faceted format: micro, cyano, geno
bray_long3 <- bray_long2 |>
  select(BrayCurtis, micro_pair, cyano_pair, geno_pair) |>
  pivot_longer(
    cols = c(micro_pair, cyano_pair, geno_pair),
    names_to = "Factor",
    values_to = "GroupType"
  )

# Rename for cleaner facet labels
bray_long3$Factor <- factor(
  bray_long3$Factor,
  levels = c("micro_pair", "cyano_pair", "geno_pair"),
  labels = c("Micro", "Cyano", "Geno")
)

# Faceted boxplot for all 3 factors
ggplot(bray_long3, aes(GroupType, BrayCurtis, fill = GroupType)) +
  geom_boxplot(outlier.size = 0.6) +
  facet_wrap(~Factor, scales = "free_x") +
  theme_bw() +
  labs(
    x = "",
    y = "Bray–Curtis Dissimilarity",
    title = "Pairwise Bray–Curtis: Within vs Between (Micro, Cyano, Geno)"
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

# MCMCglmm stats
prior <- list(
  R = list(V = 1, nu = 0.002),
  G = list(
    G1 = list(V = 1, nu = 0.002),
    G2 = list(V = 1, nu = 0.002)
  )
)

mod <- MCMCglmm(
  BrayCurtis ~ micro_pair + cyano_pair + geno_pair,
  random = ~ Sample1 + Sample2,
  data = bray_long2,
  family = "gaussian",
  prior = prior,
  nitt = 130000,
  burnin = 30000,
  thin = 100,
  verbose = FALSE
)
summary(mod)

# other stats
wilcox.test(BrayCurtis ~ micro_pair, data = bray_long2) #significant
wilcox.test(BrayCurtis ~ cyano_pair, data = bray_long2) #significant
wilcox.test(BrayCurtis ~ geno_pair, data = bray_long2) #not significant

adonis2(
  bray_dist ~ geno + cyano + micro,
  data = samples_df_filtered,
  by = "term"
)
# micro accounts for 26.3% of variance
# cyano accounts for 7.62% of variance
# geno accounts for 3.23% of variance
# 62.8% of variance is unexplained

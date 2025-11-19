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

# Cyanophyceae Analysis ------------------------------------------------------
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

# Filter out unclassified taxa before calculating relative abundance
asv_mat_long_filtered <- asv_mat_long |>
    filter(
        !is.na(Family) &
            Family != "Unclassified" &
            !is.na(Class) &
            Class != "Unclassified"
    )


# Fill missing family counts with zeroes for each treatment
asv_treatment_counts_full <- asv_mat_long_filtered |>
    complete(treatment, Family, fill = list(TotalCount = 0))


asv_treatment_counts_full$pond <- sub(
    "^([A-Za-z]+-\\d+)-.*$",
    "\\1",
    asv_treatment_counts_full$treatment
)

# Categorize treatment as 'Water' or 'Duckweed'
asv_treatment_counts_full$treatment_category <- ifelse(
    grepl("DW", asv_treatment_counts_full$treatment),
    "Duckweed",
    "Water"
)

asv_cyano <- asv_treatment_counts_full |>
    filter(Class == "Cyanophyceae")

# Diversity Analysis of Cyanophyceae -----------------------------------------

# Filter to keep only Cyanophyceae
cyano_tax_mat_filtered <- tax_mat_filtered[
    tax_mat_filtered[, "Class"] == "Cyanophyceae",
]

# Ensure unique taxa names
rownames(cyano_tax_mat_filtered) <- make.unique(rownames(
    cyano_tax_mat_filtered
))

# Filter asv matrix based on taxonomy
cyano_asv_mat_final <- asv_mat_subset_filtered[
    rownames(asv_mat_subset_filtered) %in% rownames(cyano_tax_mat_filtered),
]

# Calculate alpha diversity (Shannon index) for each sample
shannon_div <- diversity(cyano_asv_mat_final, index = "shannon", MARGIN = 2)

shannon_div_df <- data.frame(
    SampleID = colnames(cyano_asv_mat_final),
    ShannonDiversity = shannon_div
)

# Merge with sample metadata
shannon_div_df <- merge(shannon_div_df, samples_df_filtered, by = "SampleID")

# extract treatment components
shannon_div_df$pond <- sub(
    "^([A-Za-z]+-\\d+)-.*$",
    "\\1",
    shannon_div_df$treatment
)

# Categorize treatment as 'Water' or 'Duckweed'
shannon_div_df$treatment_category <- ifelse(
    grepl("DW", shannon_div_df$treatment),
    "Duckweed",
    "Water"
)
shannon_div_df$treatment_category <- factor(
    shannon_div_df$treatment_category,
    levels = c("Water", "Duckweed")
)


# Create a new pond_treatment column to control x-axis order
shannon_div_df <- shannon_div_df |>
    mutate(
        pond_treatment = paste0(
            pond,
            "_",
            ifelse(treatment_category == "Water", "W", "D")
        )
    ) |>
    mutate(
        pond_treatment = factor(
            pond_treatment,
            levels = c(
                "MP-1_W",
                "MP-1_D",
                "UM-1_W",
                "UM-1_D",
                "ODR-2_W",
                "ODR-2_D",
                "ODR-3_W",
                "ODR-3_D",
                "TF-1_W",
                "TF-1_D",
                "TF-2_W",
                "TF-2_D"
            )
        )
    ) |>
    mutate(pond_label = pond_name_mapping[pond])


# Calculate summary statistics per group
div_stats <- shannon_div_df |>
    group_by(pond, pond_label, treatment_category, pond_treatment) |>
    summarise(
        mean_div = mean(ShannonDiversity),
        sd_div = sd(ShannonDiversity),
        .groups = "drop"
    )

dodge_width <- 0.5

div_plot <- ggplot(
    div_stats,
    aes(x = pond_label, y = mean_div, color = treatment_category)
) +
    geom_jitter(
        data = shannon_div_df,
        aes(x = pond_label, y = ShannonDiversity, color = treatment_category),
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
        aes(ymin = mean_div - sd_div, ymax = mean_div + sd_div),
        width = 0.2,
        position = position_dodge(width = dodge_width)
    ) +
    scale_color_manual(
        values = c("Water" = "#88CCEE", "Duckweed" = "#44AA99"),
        name = "Sample Type"
    ) +
    labs(
        x = "Green Manure Source",
        y = "Shannon Diversity Index - Cyanophyceae"
    ) +
    scale_x_discrete(labels = unique(shannon_div_df$pond_label)) +
    theme_cowplot() +
    theme(
        legend.position = "right",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9)
    )

div_plot


#stats section

#does cyano diversity differ between water and duckweed samples
shannon_div_df |>
    group_by(treatment_category) |>
    summarise(
        mean_div = mean(ShannonDiversity),
        sd_div = sd(ShannonDiversity),
        .groups = "drop"
    )

mcmc_data <- shannon_div_df |>
    select(ShannonDiversity, pond, treatment_category) |>
    mutate(
        pond = as.factor(pond),
        treatment_category = as.factor(treatment_category)
    )

mod1 <- MCMCglmm(
    ShannonDiversity ~ treatment_category,
    data = mcmc_data,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod1) #duckweed-associated cyano are slightly significantly more diverse than water-associated cyano


#are the duckweed samples different?
mcmc_DW_data <- shannon_div_df |>
    filter(treatment_category == "Duckweed") |>
    select(ShannonDiversity, pond, treatment_category) |>
    mutate(
        pond = as.factor(pond),
        treatment_category = as.factor(treatment_category)
    )

mod2 <- MCMCglmm(
    ShannonDiversity ~ -1 + pond,
    data = mcmc_DW_data,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 1000
)
summary(mod2)

post_summ <- summary(mod2)$solutions
ci_df <- as.data.frame(post_summ)
ci_df$Effect <- rownames(ci_df)
names(ci_df)[c(1, 2, 3)] <- c("PostMean", "Lower95CI", "Upper95CI")

#slight differences in cyano on duckweeds between ponds
ggplot(ci_df, aes(x = Effect, y = PostMean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
    theme_minimal() +
    labs(
        title = "Posterior Means and 95% Credible Intervals",
        y = "Posterior Mean ± 95% CI",
        x = "Effect"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

#are the water samples different?
mcmc_W_data <- shannon_div_df |>
    filter(treatment_category == "Water") |>
    select(ShannonDiversity, pond, treatment_category) |>
    mutate(
        pond = as.factor(pond),
        treatment_category = as.factor(treatment_category)
    )

mod3 <- MCMCglmm(
    ShannonDiversity ~ -1 + pond,
    data = mcmc_W_data,
    verbose = F,
    nitt = 101000,
    thin = 10,
    burnin = 100
)
summary(mod3)

post_summ <- summary(mod3)$solutions
ci_df <- as.data.frame(post_summ)
ci_df$Effect <- rownames(ci_df)
names(ci_df)[c(1, 2, 3)] <- c("PostMean", "Lower95CI", "Upper95CI")

#slight differences in cyano in water between ponds
ggplot(ci_df, aes(x = Effect, y = PostMean)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower95CI, ymax = Upper95CI), width = 0.2) +
    theme_minimal() +
    labs(
        title = "Posterior Means and 95% Credible Intervals",
        y = "Posterior Mean ± 95% CI",
        x = "Effect"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

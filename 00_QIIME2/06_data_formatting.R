source("globals.R")

# Load environment variables
load_dot_env()

# Set paths from .env
R_inputs <- Sys.getenv("R_inputs")
output_path <- Sys.getenv("data_path")

#--------------------------------------------------
# Step 1: Prepare ASV matrix
#--------------------------------------------------

feature_table <- file.path(R_inputs, "feature_table_no_contam.tsv")

asv_mat <- read.table(
    feature_table,
    header = TRUE,
    sep = "\t",
    check.names = FALSE
)

colnames(asv_mat)[colnames(asv_mat) == "#OTU ID"] <- "asv"

asv_mat <- asv_mat |>
    tibble::column_to_rownames("asv") |>
    as.matrix()

# Write out processed ASV matrix
write.table(
    asv_mat,
    file.path(output_path, "processed_asv_matrix.tsv"),
    quote = FALSE,
    sep = "\t"
)

#--------------------------------------------------
# Step 2: Prepare taxonomy matrix
#--------------------------------------------------

tax_file <- file.path(R_inputs, "taxonomy_no_contam.tsv")

tax_mat <- read.table(
    tax_file,
    header = TRUE,
    sep = "\t"
)

colnames(tax_mat)[colnames(tax_mat) == "Feature.ID"] <- "asv"

tax_mat <- tax_mat |>
    tibble::column_to_rownames("asv") |>
    separate(
        Taxon,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        sep = ";"
    ) |>
    mutate(across(
        c(Kingdom, Phylum, Class, Order, Family, Genus),
        ~ sub("^.*?__(.*)", "\\1", .)
    ))

tax_mat <- as.matrix(tax_mat)

# Write out processed taxonomy matrix
write.table(
    tax_mat,
    file.path(output_path, "processed_taxonomy_matrix.tsv"),
    sep = "\t",
    quote = FALSE,
    col.names = NA
)

#--------------------------------------------------
# Step 3: Prepare sample data
#--------------------------------------------------

metadata_file <- file.path(R_inputs, "metadata.tsv")

samples_df <- read.table(
    metadata_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# Rename column if needed
colnames(samples_df)[colnames(samples_df) == "sampleid"] <- "sample"

# Set rownames to match ASV sample names
rownames(samples_df) <- samples_df$sample

# Now write out the processed metadata correctly
write.table(
    samples_df,
    file.path(output_path, "processed_sample_metadata.tsv"),
    quote = FALSE,
    sep = "\t",
    row.names = TRUE
)

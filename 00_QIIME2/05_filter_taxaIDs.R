source("globals.R")

# Load environment variables
load_dot_env()

# Set paths from .env
R_inputs <- Sys.getenv("R_inputs")

# ───────────────────────────────────────────────────────────────
# Load Data
# ───────────────────────────────────────────────────────────────

metadata <- read.table(
    file.path(R_inputs, "metadata.tsv"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

feature_table <- read.csv(
    file.path(R_inputs, "feature-table.csv"),
    header = TRUE,
    check.names = FALSE
)

taxonomy <- read.table(
    file.path(R_inputs, "taxonomy.tsv"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

chordata_taxa <- read.csv(
    file.path(R_inputs, "ncbi_taxid_Chordata.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

streptophyta_taxa <- read.csv(
    file.path(R_inputs, "ncbi_taxid_Streptophyta.txt"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

family_taxa_df <- read.csv(
    file.path(R_inputs, "updated_family_taxa.csv"),
    header = TRUE,
    stringsAsFactors = FALSE
)

# ───────────────────────────────────────────────────────────────
# Taxonomy Cleanup
# ───────────────────────────────────────────────────────────────
taxonomy <- taxonomy |>
    separate(
        Taxon,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        sep = ";",
        remove = TRUE
    ) |>
    mutate(across(Kingdom:Genus, ~ sub("^.*?__", "", .))) |>
    select(-Confidence)

# ───────────────────────────────────────────────────────────────
# Remove Host Contaminants (Chordata, Streptophyta)
# ───────────────────────────────────────────────────────────────
combined_taxa <- bind_rows(chordata_taxa, streptophyta_taxa)

filtered_df_family <- family_taxa_df |>
    mutate(across(c(Family, Taxid), as.character)) |>
    anti_join(
        combined_taxa |> mutate(across(c(Family, Taxid), as.character)),
        by = "Family"
    )

filtered_taxonomy_df <- taxonomy |>
    semi_join(filtered_df_family, by = "Family")

filtered_feature_table <- feature_table |>
    filter(`#OTU ID` %in% filtered_taxonomy_df$Feature.ID)

#show which ASVs were removed with their assigned taxonomy
removed_taxa <- taxonomy |>
    filter(!Feature.ID %in% filtered_taxonomy_df$Feature.ID) |>
    select(Feature.ID, Kingdom, Phylum, Class, Order, Family, Genus)

#show unique taxa removed
unique_removed_taxa <- removed_taxa |>
    group_by(Kingdom, Phylum, Class, Order, Family, Genus) |>
    summarise(Count = n()) |>
    arrange(desc(Count))

#how many were removed at this step
cat("Removed:", nrow(feature_table) - nrow(filtered_feature_table), "ASVs\n")

# ───────────────────────────────────────────────────────────────
# Identify & Remove Blanks
# ───────────────────────────────────────────────────────────────

asv_table <- filtered_feature_table |>
    column_to_rownames("#OTU ID") |>
    as.matrix() |>
    apply(2, as.numeric)

rownames(asv_table) <- filtered_feature_table$`#OTU ID`

relabund <- sweep(asv_table, 2, colSums(asv_table), `/`)

poss_blank <- rownames(relabund)[
    rowSums(relabund[, c("blank2", "blank3"), drop = FALSE]) > 0.05
]

relabund_blanks <- relabund[, grepl(
    "blank",
    colnames(relabund),
    ignore.case = TRUE
)]


#preview identified contaminants
taxonomy |> filter(Feature.ID %in% poss_blank)
round(relabund[poss_blank, ], 2)

#how many were removed at this step
cat("Removed:", length(poss_blank), "OTUs\n")

# ───────────────────────────────────────────────────────────────
# Final Cleanup: Remove Contaminants
# ───────────────────────────────────────────────────────────────
feature_table_no_contam <- filtered_feature_table |>
    filter(!`#OTU ID` %in% poss_blank)

taxonomy_no_contam <- filtered_taxonomy_df |>
    filter(!Feature.ID %in% poss_blank)

# Sanity check
cat("Removed:", length(poss_blank), "contaminants\n")


# ───────────────────────────────────────────────────────────────
# Save Cleaned Files
# ───────────────────────────────────────────────────────────────

# Make sure first column is named consistently
colnames(feature_table_no_contam)[1] <- "OTU_ID"

# Build output file paths within R_inputs
feature_out <- file.path(R_inputs, "feature_table_no_contam.tsv")
taxonomy_out <- file.path(R_inputs, "taxonomy_no_contam.tsv")

# Write out cleaned feature table
write.table(
    feature_table_no_contam,
    feature_out,
    sep = "\t",
    row.names = TRUE,
    quote = FALSE
)

# Write out cleaned taxonomy
write.table(
    taxonomy_no_contam,
    taxonomy_out,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

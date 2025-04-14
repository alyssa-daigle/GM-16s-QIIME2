library(tidyr)
library(dplyr)
library(decontam)
library(tibble)

# ───────────────────────────────────────────────────────────────
# 1. Setup
# ───────────────────────────────────────────────────────────────
setwd("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/16s_Analysis/taxize")

# ───────────────────────────────────────────────────────────────
# 2. Load Data
# ───────────────────────────────────────────────────────────────
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
feature_table <- read.csv("feature-table.csv", header = TRUE, check.names = FALSE)
taxonomy <- read.table("taxonomy.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chordata_taxa <- read.csv("ncbi_taxid_Chordata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
streptophyta_taxa <- read.csv("ncbi_taxid_Streptophyta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
family_taxa_df <- read.csv("family_taxa.csv", header = TRUE, stringsAsFactors = FALSE)

# ───────────────────────────────────────────────────────────────
# 3. Taxonomy Cleanup
# ───────────────────────────────────────────────────────────────
taxonomy <- taxonomy |> 
    separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = TRUE) |> 
    mutate(across(Kingdom:Genus, ~sub("^.*?__", "", .))) |> 
    select(-Confidence)

# ───────────────────────────────────────────────────────────────
# 4. Remove Host Contaminants (Chordata, Streptophyta)
# ───────────────────────────────────────────────────────────────
combined_taxa <- bind_rows(chordata_taxa, streptophyta_taxa)

filtered_df_family <- family_taxa_df |> 
    mutate(across(c(Family, Taxid), as.character)) |> 
    anti_join(combined_taxa |>  mutate(across(c(Family, Taxid), as.character)), by = "Family")

filtered_taxonomy_df <- taxonomy |> 
    semi_join(filtered_df_family, by = "Family")

filtered_feature_table <- feature_table |> 
    filter(`#OTU ID` %in% filtered_taxonomy_df$Feature.ID)

# ───────────────────────────────────────────────────────────────
# 5. Identify & Remove Blanks 
# ───────────────────────────────────────────────────────────────

otu_table <- filtered_feature_table |> 
    column_to_rownames("#OTU ID") |> 
    as.matrix() |> 
    apply(2, as.numeric)

rownames(otu_table) <- filtered_feature_table$`#OTU ID`

relabund <- sweep(otu_table, 2, colSums(otu_table), `/`)

poss_blank <- rownames(relabund)[
    rowSums(relabund[, c("blank2", "blank3"), drop = FALSE]) > 0.05
]

# Preview identified contaminants
taxonomy |>  filter(Feature.ID %in% poss_blank)
round(relabund[poss_blank, ], 1)

# ───────────────────────────────────────────────────────────────
# 6. Final Cleanup: Remove Contaminants
# ───────────────────────────────────────────────────────────────
feature_table_no_contam <- filtered_feature_table |> 
    filter(!`#OTU ID` %in% poss_blank)

taxonomy_no_contam <- filtered_taxonomy_df |> 
    filter(!Feature.ID %in% poss_blank)

# Sanity check
cat("Removed:", length(poss_blank), "contaminants\n")

# ───────────────────────────────────────────────────────────────
# 7. Save Cleaned Files
# ───────────────────────────────────────────────────────────────
colnames(feature_table_no_contam)[1] <- "OTU_ID"
write.table(feature_table_no_contam, "feature_table_no_contam.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

write.table(taxonomy_no_contam, "taxonomy_no_contam.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

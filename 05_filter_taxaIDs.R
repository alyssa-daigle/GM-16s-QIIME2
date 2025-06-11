library(tidyr)
library(dplyr)
library(tibble)
library(MCMCglmm)
library(data.table)
library(ggtext) 
library(cowplot)
library(ggplot2)

# ───────────────────────────────────────────────────────────────
# Load Data
# ───────────────────────────────────────────────────────────────
metadata <- read.table("metadata.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
feature_table <- read.csv("feature-table.csv", header = TRUE, check.names = FALSE)
taxonomy <- read.table("taxonomy.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chordata_taxa <- read.csv("ncbi_taxid_Chordata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
streptophyta_taxa <- read.csv("ncbi_taxid_Streptophyta.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
family_taxa_df <- read.csv("updated_family_taxa.csv", header = TRUE, stringsAsFactors = FALSE)

# ───────────────────────────────────────────────────────────────
# Taxonomy Cleanup
# ───────────────────────────────────────────────────────────────
taxonomy <- taxonomy |> 
    separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";", remove = TRUE) |> 
    mutate(across(Kingdom:Genus, ~sub("^.*?__", "", .))) |> 
    select(-Confidence)

# ───────────────────────────────────────────────────────────────
# Remove Host Contaminants (Chordata, Streptophyta)
# ───────────────────────────────────────────────────────────────
combined_taxa <- bind_rows(chordata_taxa, streptophyta_taxa)

filtered_df_family <- family_taxa_df |> 
    mutate(across(c(Family, Taxid), as.character)) |> 
    anti_join(combined_taxa |>  mutate(across(c(Family, Taxid), as.character)), by = "Family")

filtered_taxonomy_df <- taxonomy |> 
    semi_join(filtered_df_family, by = "Family")

filtered_feature_table <- feature_table |> 
    filter(`#OTU ID` %in% filtered_taxonomy_df$Feature.ID)

#show which OTUs were removed with their assigned taxonomy
removed_taxa <- taxonomy |> 
    filter(!Feature.ID %in% filtered_taxonomy_df$Feature.ID) |> 
    select(Feature.ID, Kingdom, Phylum, Class, Order, Family, Genus)

#show unique taxa removed
unique_removed_taxa <- removed_taxa |> 
    group_by(Kingdom, Phylum, Class, Order, Family, Genus) |> 
    summarise(Count = n()) |> 
    arrange(desc(Count))

#how many were removed at this step
cat("Removed:", nrow(feature_table) - nrow(filtered_feature_table), "OTUs\n")

# ───────────────────────────────────────────────────────────────
# Identify & Remove Blanks 
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

relabund_blanks <- relabund[, grepl("blank", colnames(relabund), ignore.case = TRUE)]


#preview identified contaminants
taxonomy |>  filter(Feature.ID %in% poss_blank)
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
colnames(feature_table_no_contam)[1] <- "OTU_ID"
write.table(feature_table_no_contam, "feature_table_no_contam.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

write.table(taxonomy_no_contam, "taxonomy_no_contam.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# ───────────────────────────────────────────────────────────────
# Summary Statistics (Expt 1)
# ───────────────────────────────────────────────────────────────

#define the subset of samples to include
selected_samples <- colnames(feature_table)[-1]  # remove OTU ID
selected_samples <- selected_samples[selected_samples <= "24-d6"]

#filter columns for selected samples (before filtering)
feature_table_subset <- feature_table[, c("#OTU ID", selected_samples)]

#filter columns for selected samples (after filtering)
feature_table_no_contam_subset <- feature_table_no_contam[, c("OTU_ID", selected_samples)]

#calculate OTU counts
total_otus_start <- nrow(feature_table_subset)
total_otus_removed <- total_otus_start - nrow(feature_table_no_contam_subset)
total_otus_remaining <- nrow(feature_table_no_contam_subset)

#calculate total reads
total_reads_start <- sum(feature_table_subset[, -1])
total_reads_end <- sum(feature_table_no_contam_subset[, -1])
total_reads_removed <- total_reads_start - total_reads_end

#per-sample read stats (before filtering)
reads_per_sample_start <- colSums(feature_table_subset[, -1])
max_reads_per_sample_start <- max(reads_per_sample_start)
min_reads_per_sample_start <- min(reads_per_sample_start)
median_reads_per_sample_start <- median(reads_per_sample_start)
mean_reads_per_sample_start <- mean(reads_per_sample_start)

#per-sample read stats (after filtering)
reads_per_sample_end <- colSums(feature_table_no_contam_subset[, -1])
max_reads_per_sample_end <- max(reads_per_sample_end)
min_reads_per_sample_end <- min(reads_per_sample_end)
median_reads_per_sample_end <- median(reads_per_sample_end)
mean_reads_per_sample_end <- mean(reads_per_sample_end)

cat(
    "Summary for Expt 1 sample:\n",
    "--------------------------------------------------\n",
    "Total OTUs at start:", total_otus_start, "\n",
    "Total OTUs removed:", total_otus_removed, "\n",
    "Total OTUs remaining:", total_otus_remaining, "\n",
    "Total reads at start:", total_reads_start, "\n",
    "Total reads at end:", total_reads_end, "\n",
    "Total reads removed:", total_reads_removed, "\n\n",
    "Reads per sample BEFORE filtering:\n",
    "Max:", max_reads_per_sample_start, "\n",
    "Min:", min_reads_per_sample_start, "\n",
    "Median:", median_reads_per_sample_start, "\n",
    "Mean:", mean_reads_per_sample_start, "\n\n",
    "Reads per sample AFTER filtering:\n",
    "Max:", max_reads_per_sample_end, "\n",
    "Min:", min_reads_per_sample_end, "\n",
    "Median:", median_reads_per_sample_end, "\n",
    "Mean:", mean_reads_per_sample_end, "\n"
)


# ───────────────────────────────────────────────────────────────
# Summary Statistics (Expt 2)
# ───────────────────────────────────────────────────────────────
# Get all sample names from the original filtered_feature_table (before contaminant removal)
all_samples <- colnames(filtered_feature_table)[-1]  # exclude OTU ID column
samples_excluded <- all_samples[1:192]
selected_samples<- all_samples[!(all_samples %in% samples_excluded |
                                     grepl("blank", all_samples, ignore.case = TRUE) |
                                     grepl("Expt1$", all_samples) |
                                     grepl("Expt2$", all_samples))]

#filter columns for selected samples (before filtering)
feature_table_subset <- feature_table[, c("#OTU ID", selected_samples)]

#filter columns for selected samples (after filtering)
feature_table_no_contam_subset <- feature_table_no_contam[, c("OTU_ID", selected_samples)]

#calculate OTU counts
total_otus_start <- nrow(feature_table_subset)
total_otus_removed <- total_otus_start - nrow(feature_table_no_contam_subset)
total_otus_remaining <- nrow(feature_table_no_contam_subset)

#calculate total reads
total_reads_start <- sum(feature_table_subset[, -1])
total_reads_end <- sum(feature_table_no_contam_subset[, -1])
total_reads_removed <- total_reads_start - total_reads_end

#per-sample read stats (before filtering)
reads_per_sample_start <- colSums(feature_table_subset[, -1])
max_reads_per_sample_start <- max(reads_per_sample_start)
min_reads_per_sample_start <- min(reads_per_sample_start)
median_reads_per_sample_start <- median(reads_per_sample_start)
mean_reads_per_sample_start <- mean(reads_per_sample_start)

#per-sample read stats (after filtering)
reads_per_sample_end <- colSums(feature_table_no_contam_subset[, -1])
max_reads_per_sample_end <- max(reads_per_sample_end)
min_reads_per_sample_end <- min(reads_per_sample_end)
median_reads_per_sample_end <- median(reads_per_sample_end)
mean_reads_per_sample_end <- mean(reads_per_sample_end)

cat(
    "Summary for Expt 2 Samples:\n",
    "--------------------------------------------------\n",
    "Total OTUs at start:", total_otus_start, "\n",
    "Total OTUs removed:", total_otus_removed, "\n",
    "Total OTUs remaining:", total_otus_remaining, "\n",
    "Total reads at start:", total_reads_start, "\n",
    "Total reads at end:", total_reads_end, "\n",
    "Total reads removed:", total_reads_removed, "\n\n",
    "Reads per sample BEFORE filtering:\n",
    "Max:", max_reads_per_sample_start, "\n",
    "Min:", min_reads_per_sample_start, "\n",
    "Median:", median_reads_per_sample_start, "\n",
    "Mean:", mean_reads_per_sample_start, "\n\n",
    "Reads per sample AFTER filtering:\n",
    "Max:", max_reads_per_sample_end, "\n",
    "Min:", min_reads_per_sample_end, "\n",
    "Median:", median_reads_per_sample_end, "\n",
    "Mean:", mean_reads_per_sample_end, "\n"
)


# ───────────────────────────────────────────────────────────────
# Test for differences in read counts across micro treatments
# ───────────────────────────────────────────────────────────────

col_index <- which(colnames(feature_table_no_contam) == "24-d6")
feature_table_no_contam <- feature_table_no_contam[, 1:col_index]

rownames(feature_table_no_contam) <- feature_table_no_contam[, 1]
feature_table_no_contam <- feature_table_no_contam[, -1]

#pivot the table longer
cleaned_feature_table_long <- feature_table_no_contam %>%
    rownames_to_column("feature_id") %>%
    pivot_longer(cols = -feature_id, names_to = "sampleid", values_to = "count")

clean_merged_table <- cleaned_feature_table_long %>%
    left_join(metadata, by = "sampleid")

#convert to data.table for faster processing
setDT(clean_merged_table)

#split the 'treatment' column and assign to the appropriate columns
clean_merged_table[, c("geno", "cyano", "micro") := tstrsplit(treatment, "-", fixed = TRUE)]

#ceate inoculation_status using vectorized logic
clean_merged_table[, inoculation_status := fcase(
    cyano == "N" & micro == "N", "Uninoculated",
    cyano == "Y" & micro == "N", "*M. aeruginosa*<br> only",
    cyano == "N" & micro %in% c("H", "ODR", "KF"), "Microbiome<br> only",
    cyano == "Y" & micro %in% c("H", "ODR", "KF"), "*M. aeruginosa*<br>+ Microbiome"
)]

clean_merged_table[, inoculation_status := factor(inoculation_status,
                                                  levels = c("Uninoculated", "*M. aeruginosa*<br> only", "Microbiome<br> only", "*M. aeruginosa*<br>+ Microbiome"
)
)]

#remove 0 counts
clean_merged_table <- clean_merged_table %>%
    filter(count > 0)

#make a table to summarize number of reads in inoculated, uninoculated, and other to compare
clean_merged_table_summary <- clean_merged_table %>%
    group_by(inoculation_status) %>%
    summarise(total_reads = sum(count, na.rm = TRUE),
              num_samples = n()) %>%
    arrange(desc(total_reads))

#summarize total reads per sample
sample_level <- clean_merged_table[, .(sample_reads = sum(count)), by = .(sampleid, inoculation_status)]

#mean and SD per group
summary_table <- sample_level[, .(
    mean_reads_per_sample = mean(sample_reads),
    sd_reads_per_sample = sd(sample_reads),
    n_samples = .N
), by = inoculation_status]


#model to determine significant differences
model <- MCMCglmm(sample_reads ~ -1 + inoculation_status,
                  data = sample_level,
                  nitt = 101000, burnin = 1000, thin = 10,
                  verbose = FALSE)
summary(model)

# Create a new column for significance letters
summary_table$significance <- c("c", "bc", "a", "ab")

#make a list of colors for each inoculation_status in summary_table
# Get unique inoculation statuses
statuses <- unique(summary_table$inoculation_status)

# Now, plot with annotation
ggplot(summary_table, aes(x = inoculation_status, y = mean_reads_per_sample)) +
    geom_col(show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean_reads_per_sample - sd_reads_per_sample,
                      ymax = mean_reads_per_sample + sd_reads_per_sample),
                  width = 0.2) +
    geom_text(aes(y = mean_reads_per_sample + sd_reads_per_sample, label = significance), 
              vjust = -0.5, size = 5, color = "black") +
    labs(
        x = "Inoculation Status",
        y = "Mean Reads Per Sample"
    ) +
    theme_cowplot() +
    scale_y_continuous(labels = scales::comma, 
                       breaks = seq(0, max(summary_table$mean_reads_per_sample + summary_table$sd_reads_per_sample), by = 25000), 
                       expand = expansion(mult = c(0, 0.05))) +
    theme(axis.text.x = element_markdown(size = 9))

ggsave("Expt1_mean_reads_per_sample.jpg", width = 7, height = 6, dpi = 500)
ggsave("~/Library/CloudStorage/OneDrive-UniversityofNewHampshire/GreenManureProject/WRITING/plots/Expt1_mean_reads_per_sample.jpg", width = 7, height = 6, dpi = 500)


## metadata --------------------------------------------
# Split the 'treatment' column and assign to the appropriate columns

#subset metadata up to row 192
metadata <- metadata[1:192, ]

metadata <- metadata |> 
  separate(treatment, into = c("geno", "cyano", "micro"), sep = "-", remove = FALSE) |> 
  mutate(
    inoculation_status = case_when(
      cyano == "N" & micro == "N" ~ "Uninoculated",
      cyano == "Y" & micro == "N" ~ "*M. aeruginosa*<br> only",
      cyano == "N" & micro %in% c("H", "ODR", "KF") ~ "Microbiome<br> only",
      cyano == "Y" & micro %in% c("H", "ODR", "KF") ~ "*M. aeruginosa*<br>+ Microbiome"
    )
  )

#count number of samples per inoculation status
metadata_summary <- metadata |> 
  group_by(inoculation_status) |> 
  summarise(num_samples = n()) |> 
  arrange(desc(num_samples))

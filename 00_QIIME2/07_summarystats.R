# Load environment variables
source("globals.R")
load_dot_env()

# Set paths from .env
R_inputs <- Sys.getenv("R_inputs")

# ───────────────────────────────────────────────────────────────
# Load required data
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

feature_table_no_contam <- read.table(
  file.path(R_inputs, "feature_table_no_contam.tsv"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

filtered_feature_table <- read.csv(
  file.path(R_inputs, "filtered_feature_table.csv"),
  header = TRUE,
  check.names = FALSE
)

# ───────────────────────────────────────────────────────────────
# Summary Statistics (Expt 1)
# ───────────────────────────────────────────────────────────────

# define the subset of samples to include
selected_samples <- colnames(feature_table)[-1] # remove ASV ID
selected_samples <- selected_samples[selected_samples <= "24-d6"]

# filter columns for selected samples (before filtering)
feature_table_subset <- feature_table[, c("#OTU ID", selected_samples)]

# filter columns for selected samples (after filtering)
feature_table_no_contam_subset <- feature_table_no_contam[, c(
  "OTU_ID",
  selected_samples
)]

# calculate ASV counts
total_asvs_start <- nrow(feature_table_subset)
total_asvs_removed <- total_asvs_start - nrow(feature_table_no_contam_subset)
total_asvs_remaining <- nrow(feature_table_no_contam_subset)

# calculate total reads
total_reads_start <- sum(feature_table_subset[, -1])
total_reads_end <- sum(feature_table_no_contam_subset[, -1])
total_reads_removed <- total_reads_start - total_reads_end

# per-sample read stats (before filtering)
reads_per_sample_start <- colSums(feature_table_subset[, -1])
max_reads_per_sample_start <- max(reads_per_sample_start)
min_reads_per_sample_start <- min(reads_per_sample_start)
median_reads_per_sample_start <- median(reads_per_sample_start)
mean_reads_per_sample_start <- mean(reads_per_sample_start)

# per-sample read stats (after filtering)
reads_per_sample_end <- colSums(feature_table_no_contam_subset[, -1])
max_reads_per_sample_end <- max(reads_per_sample_end)
min_reads_per_sample_end <- min(reads_per_sample_end)
median_reads_per_sample_end <- median(reads_per_sample_end)
mean_reads_per_sample_end <- mean(reads_per_sample_end)

cat(
  "Summary for Expt 1 sample:\n",
  "--------------------------------------------------\n",
  "Total ASVs at start:",
  total_asvs_start,
  "\n",
  "Total ASVs removed:",
  total_asvs_removed,
  "\n",
  "Total ASVs remaining:",
  total_asvs_remaining,
  "\n",
  "Total reads at start:",
  total_reads_start,
  "\n",
  "Total reads at end:",
  total_reads_end,
  "\n",
  "Total reads removed:",
  total_reads_removed,
  "\n\n",
  "Reads per sample BEFORE filtering:\n",
  "Max:",
  max_reads_per_sample_start,
  "\n",
  "Min:",
  min_reads_per_sample_start,
  "\n",
  "Median:",
  median_reads_per_sample_start,
  "\n",
  "Mean:",
  mean_reads_per_sample_start,
  "\n\n",
  "Reads per sample AFTER filtering:\n",
  "Max:",
  max_reads_per_sample_end,
  "\n",
  "Min:",
  min_reads_per_sample_end,
  "\n",
  "Median:",
  median_reads_per_sample_end,
  "\n",
  "Mean:",
  mean_reads_per_sample_end,
  "\n"
)

# ───────────────────────────────────────────────────────────────
# Summary Statistics (Expt 2)
# ───────────────────────────────────────────────────────────────

all_samples <- colnames(filtered_feature_table)[-1] # exclude ASV ID
samples_excluded <- all_samples[1:192]

selected_samples <- all_samples[
  !(all_samples %in%
    samples_excluded |
    grepl("blank", all_samples, ignore.case = TRUE) |
    grepl("Expt1$", all_samples) |
    grepl("Expt2$", all_samples))
]

feature_table_subset <- feature_table[, c("#OTU ID", selected_samples)]
feature_table_no_contam_subset <- feature_table_no_contam[, c(
  "OTU_ID",
  selected_samples
)]

total_asvs_start <- nrow(feature_table_subset)
total_asvs_removed <- total_asvs_start - nrow(feature_table_no_contam_subset)
total_asvs_remaining <- nrow(feature_table_no_contam_subset)

total_reads_start <- sum(feature_table_subset[, -1])
total_reads_end <- sum(feature_table_no_contam_subset[, -1])
total_reads_removed <- total_reads_start - total_reads_end

reads_per_sample_start <- colSums(feature_table_subset[, -1])
reads_per_sample_end <- colSums(feature_table_no_contam_subset[, -1])

cat(
  "Summary for Expt 2 Samples:\n",
  "--------------------------------------------------\n",
  "Total ASVs at start:",
  total_asvs_start,
  "\n",
  "Total ASVs removed:",
  total_asvs_removed,
  "\n",
  "Total ASVs remaining:",
  total_asvs_remaining,
  "\n",
  "Total reads at start:",
  total_reads_start,
  "\n",
  "Total reads at end:",
  total_reads_end,
  "\n",
  "Total reads removed:",
  total_reads_removed,
  "\n"
)

# ───────────────────────────────────────────────────────────────
# Test for differences in read counts across micro treatments
# ───────────────────────────────────────────────────────────────

col_index <- which(colnames(feature_table_no_contam) == "24-d6")
feature_table_no_contam <- feature_table_no_contam[, 1:col_index]

rownames(feature_table_no_contam) <- feature_table_no_contam[, 1]
feature_table_no_contam <- feature_table_no_contam[, -1]

#pivot the table longer
cleaned_feature_table_long <- feature_table_no_contam |>
  rownames_to_column("feature_id") |>
  pivot_longer(cols = -feature_id, names_to = "sampleid", values_to = "count")

clean_merged_table <- cleaned_feature_table_long |>
  left_join(metadata, by = "sampleid")

#convert to data.table for faster processing
setDT(clean_merged_table)

#split the 'treatment' column and assign to the appropriate columns
clean_merged_table[,
  c("geno", "cyano", "micro") := tstrsplit(treatment, "-", fixed = TRUE)
]

#ceate inoculation_status using vectorized logic
clean_merged_table[,
  inoculation_status := fcase(
    cyano == "N" & micro == "N"                   , "Uninoculated"                    ,
    cyano == "Y" & micro == "N"                   , "*M. aeruginosa*<br> only"        ,
    cyano == "N" & micro %in% c("H", "ODR", "KF") , "Microbiome<br> only"             ,
    cyano == "Y" & micro %in% c("H", "ODR", "KF") , "*M. aeruginosa*<br>+ Microbiome"
  )
]

clean_merged_table[,
  inoculation_status := factor(
    inoculation_status,
    levels = c(
      "Uninoculated",
      "*M. aeruginosa*<br> only",
      "Microbiome<br> only",
      "*M. aeruginosa*<br>+ Microbiome"
    )
  )
]

#remove 0 counts
clean_merged_table <- clean_merged_table |>
  filter(count > 0)

#make a table to summarize number of reads in inoculated, uninoculated, and other to compare
clean_merged_table_summary <- clean_merged_table |>
  group_by(inoculation_status) |>
  summarise(total_reads = sum(count, na.rm = TRUE), num_samples = n()) |>
  arrange(desc(total_reads))

#summarize total reads per sample
sample_level <- clean_merged_table[,
  .(sample_reads = sum(count)),
  by = .(sampleid, inoculation_status)
]

#mean and SD per group
summary_table <- sample_level[,
  .(
    mean_reads_per_sample = mean(sample_reads),
    sd_reads_per_sample = sd(sample_reads),
    n_samples = .N
  ),
  by = inoculation_status
]


#model to determine significant differences
model <- MCMCglmm(
  sample_reads ~ -1 + inoculation_status,
  data = sample_level,
  nitt = 101000,
  burnin = 1000,
  thin = 10,
  verbose = FALSE
)
summary(model)

# Create a new column for significance letters
summary_table$significance <- c("c", "bc", "a", "ab")

#make a list of colors for each inoculation_status in summary_table
# Get unique inoculation statuses
statuses <- unique(summary_table$inoculation_status)

# Now, plot with annotation
ggplot(summary_table, aes(x = inoculation_status, y = mean_reads_per_sample)) +
  geom_col(show.legend = FALSE) +
  geom_errorbar(
    aes(
      ymin = mean_reads_per_sample - sd_reads_per_sample,
      ymax = mean_reads_per_sample + sd_reads_per_sample
    ),
    width = 0.2
  ) +
  geom_text(
    aes(
      y = mean_reads_per_sample + sd_reads_per_sample,
      label = significance
    ),
    vjust = -0.5,
    size = 5,
    color = "black"
  ) +
  labs(
    x = "Inoculation Status",
    y = "Mean Reads Per Sample"
  ) +
  theme_cowplot() +
  scale_y_continuous(
    labels = scales::comma,
    breaks = seq(
      0,
      max(
        summary_table$mean_reads_per_sample +
          summary_table$sd_reads_per_sample
      ),
      by = 25000
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme(axis.text.x = element_markdown(size = 9))

ggsave("Expt1_mean_reads_per_sample.jpg", width = 7, height = 6, dpi = 500)


## metadata --------------------------------------------
# Split the 'treatment' column and assign to the appropriate columns

#subset metadata up to row 192
metadata <- metadata[1:192, ]

metadata <- metadata |>
  separate(
    treatment,
    into = c("geno", "cyano", "micro"),
    sep = "-",
    remove = FALSE
  ) |>
  mutate(
    inoculation_status = case_when(
      cyano == "N" & micro == "N" ~ "Uninoculated",
      cyano == "Y" & micro == "N" ~ "*M. aeruginosa*<br> only",
      cyano == "N" &
        micro %in% c("H", "ODR", "KF") ~ "Microbiome<br> only",
      cyano == "Y" &
        micro %in%
          c("H", "ODR", "KF") ~ "*M. aeruginosa*<br>+ Microbiome"
    )
  )

#count number of samples per inoculation status
metadata_summary <- metadata |>
  group_by(inoculation_status) |>
  summarise(num_samples = n()) |>
  arrange(desc(num_samples))

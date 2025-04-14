library(tidyr)
library(taxize)
library(purrr)
library(dplyr)

#loading taxonomy data
taxonomy <- read.table("~/Desktop/taxize/taxonomy.tsv", 
                       header = TRUE, 
                       sep = "\t", 
                       stringsAsFactors = FALSE)


taxonomy <- taxonomy |> 
    separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), 
             sep = ";", 
             remove = FALSE)

taxonomy$Kingdom <- sub("^.*?__(.*)", "\\1", taxonomy$Kingdom)
taxonomy$Phylum <- sub("^.*?__(.*)", "\\1", taxonomy$Phylum)
taxonomy$Class <- sub("^.*?__(.*)", "\\1", taxonomy$Class)
taxonomy$Order <- sub("^.*?__(.*)", "\\1", taxonomy$Order)
taxonomy$Family <- sub("^.*?__(.*)", "\\1", taxonomy$Family)
taxonomy$Genus <- sub("^.*?__(.*)", "\\1", taxonomy$Genus)

#filtering by family
unique_families <- unique(taxonomy$Family)

taxize::taxize_options(ncbi_sleep = 0.5)

filter = "^(bacteria.*|.*phyta$|.*archaeota$)"
family_ncbi_ids <- get_uid(unique_families, 
                           rank_query = "family",
                           division_filter = filter)


family_taxa_df <- tibble(Family = unique_families, Taxid = family_ncbi_ids)


missing_family_df <- family_taxa_df |> 
	filter(is.na(Taxid))


#write out as csv
write.csv(family_taxa_df, "~/Desktop/family_taxa.csv", row.names = FALSE)





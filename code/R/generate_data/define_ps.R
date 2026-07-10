# Clear environment
rm(list = ls())

library(qiime2R)
library(phyloseq)
library(tidyverse)
source("./code/R/generate_data/ps_agglom_function.R")

# define sample names
size <- data.frame(
  ranges = c("0-0.85","0.85-1.4", "1.4-2", "2-2.8", "2.8-4", ">4"),
  name = c("Floccular", "S", "M", "L", "XL", "XXL"),
  midpoint = c(0.425, 1.125, 1.7, 2.4, 3.4, 4.5)
)

# Import QIIME2 data as phyloseq object
ps <- qiime2R::qza_to_phyloseq(
  features = "./data/qiime/table_dada2.qza",
  tree = "./data/qiime/rooted_tree.qza",
  taxonomy = "./data/qiime/taxonomy.qza",
  metadata = "./data/qiime/sample-metadata.tsv"
)

ps@sam_data$size.mm       <- factor(ps@sam_data$size.mm, levels = size$ranges)
ps@sam_data$size.name     <- factor(size$name[as.numeric(ps@sam_data$size.mm)], levels = size$name)
ps@sam_data$size.midpoint <- size$midpoint[as.numeric(ps@sam_data$size.mm)]
  
# ------ Filter ------

remove_names <- taxa_names(
  subset_taxa(
    ps,
    Kingdom == "Unassigned" |    
      Order == "Chloroplast" |
      Family == "Mitochondria"
  )
)

keep_taxa <- !(taxa_names(ps) %in% remove_names)

ps_filt <- prune_taxa(keep_taxa, ps)

# ------ Rarefy ------

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps_filt))  # lowest number of ASVs per sample

ps_rare <- rarefy_even_depth(
  ps_filt, rarefy_level, rngseed = 1, replace = FALSE, trimOTUs = TRUE, verbose = TRUE
)

# ------ Save at ASV level ------

saveRDS(ps_filt, file = "./data/phyloseq/ps_ASV.rds")
saveRDS(ps_rare, file = "./data/phyloseq/ps_ASV_rarefied.rds")

# ------ Agglomerate, keeping NA values  ------

ps_genus   <- agglom_genus(ps_filt)
ps_species <- agglom_species(ps_filt)

# ------ Save at genus level ------

# remove floccular granules
ps_sub <- subset_samples(ps_genus, size.name != "Floccular")

saveRDS(ps_genus, file = "./data/phyloseq/ps_genus_full.rds")
saveRDS(ps_sub,   file = "./data/phyloseq/ps_genus_subset.rds")

# ------ Species level ------

saveRDS(ps_species, file = "./data/phyloseq/ps_species_full.rds")
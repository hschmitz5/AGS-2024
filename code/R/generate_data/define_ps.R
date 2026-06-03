# Clear environment
rm(list = ls())

library(qiime2R)
library(phyloseq)
library(tidyverse)

# define sample names
size <- data.frame(
  ranges = c("0.43-0.85","0.85-1.4", "1.4-2", "2-2.8", "2.8-4", ">4"),
  name = c("floccular", "S", "M", "L", "XL", "XXL")
)

# Import QIIME2 data as phyloseq object
ps <- qiime2R::qza_to_phyloseq(
  features = "./data/qiime/table_dada2.qza",
  tree = "./data/qiime/rooted_tree.qza",
  taxonomy = "./data/qiime/taxonomy.qza",
  metadata = "./data/qiime/sample-metadata.tsv"
)

ps@sam_data$size.mm <- factor(ps@sam_data$size.mm, levels = size$ranges)
ps@sam_data$size.name <- factor(size$name[as.numeric(ps@sam_data$size.mm)], levels = size$name)
  
# ------ Filter ------

# remove Mitochondria and Chloroplasts (removes Eukaryotes)
ps_filt0 <- phyloseq::subset_taxa(ps, ! Family %in% c("Mitochondria", "Chloroplast"))
# remove unclassified sequences
ps_filt0 <- phyloseq::subset_taxa(ps, Kingdom != "Unassigned")
  
# ------ Rarefy ------

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps_filt0))  # lowest number of ASVs per sample
# apply rarefaction
ps_filt <-rarefy_even_depth(
  ps_filt0, rarefy_level, rngseed = 7, replace = TRUE, trimOTUs = TRUE, verbose = TRUE
)

# ------ Generate new names for OTUs -------

taxonomy <- as.data.frame(as.matrix(ps_filt@tax_table)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Species_tmp = case_when(
      is.na(Phylum) ~ paste0("Phylum_unknown"),
      !is.na(Species) & !startsWith(Species, "midas") ~ paste0(Species, "_s"), 
      !is.na(Genus)   & !startsWith(Genus, "midas")   ~ paste0(Genus,   "_g"),
      !is.na(Family)  & !startsWith(Family, "midas")  ~ paste0(Family,  "_f"),
      !is.na(Order)   & !startsWith(Order, "midas")   ~ paste0(Order,   "_o"),
      !is.na(Class)   & !startsWith(Class, "midas")   ~ paste0(Class,   "_c"),
      !is.na(Phylum)  & !startsWith(Phylum, "midas")  ~ paste0(Phylum,  "_p"),
      .default = Species 
    )
  ) %>%
  group_by(Species_tmp) %>%
  mutate(
    Species_updated = if (n() == 1) {
      Species_tmp
    } else {
      paste0(Species_tmp, "-", row_number())
    }
  ) %>%
  ungroup() %>%
  dplyr::select(-Species_tmp) %>%
  column_to_rownames("OTU") %>%
  as.matrix()

# Add Species_updated to phyloseq object
tax_table(ps_filt) <- tax_table(taxonomy)
  
# Save
saveRDS(ps_filt, file = "./data/ps_ASV_full.rds")

# remove XS granules
ps_sub <- subset_samples(ps_filt, size.mm != "0.43-0.85")
  
saveRDS(ps_sub, file = "./data/ps_ASV_subset.rds")

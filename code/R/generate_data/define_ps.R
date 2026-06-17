# Clear environment
rm(list = ls())

library(qiime2R)
library(phyloseq)
library(tidyverse)

# define sample names
size <- data.frame(
  ranges = c("0-0.85","0.85-1.4", "1.4-2", "2-2.8", "2.8-4", ">4"),
  name = c("floccular", "S", "M", "L", "XL", "XXL"),
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

ps_filt <- phyloseq::subset_taxa(ps,
                                 # remove unclassified sequences
                                 (Kingdom != "Unassigned") &
                                 # remove Chloroplasts and Mitochondria
                                 (Order != "Chloroplast") &
                                 (Family != "Mitochondria")
)

# ------ Rarefy ------

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps_filt))  # lowest number of ASVs per sample

ps_rarefied <- rarefy_even_depth(
  ps_filt, rarefy_level, rngseed = 1, replace = FALSE, trimOTUs = TRUE, verbose = TRUE
)

# ------ Save at ASV level ------

saveRDS(ps_filt,     file = "./data/phyloseq/ps_ASV.rds")
saveRDS(ps_rarefied, file = "./data/phyloseq/ps_ASV_rarefied.rds")

# ------ Agglomerate, keeping NA values  ------

# If Genus is NA, then replace with higher order
taxonomy <- data.frame(ps_filt@tax_table) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Genus = case_when(
      !is.na(Genus)  ~ Genus,
      !is.na(Family) ~ ifelse(startsWith(Family, "midas"), paste0("Unk_", Family), paste0("Unk_f_", Family)),
      !is.na(Order)  ~ ifelse(startsWith(Order, "midas"), paste0("Unk_", Order), paste0("Unk_o_", Order)),
      !is.na(Class)  ~ ifelse(startsWith(Class, "midas"), paste0("Unk_", Class), paste0("Unk_c_", Class)),
      !is.na(Phylum) ~ ifelse(startsWith(Phylum, "midas"), paste0("Unk_", Phylum), paste0("Unk_p_", Phylum)),
      is.na(Phylum)  ~ paste0("Unknown_Phylum"),
      .default = Genus
    )
  ) %>%
  # name ASVs
  group_by(Genus) %>%
  mutate(
    new_OTU = if (n() == 1) {
      Genus
    } else {
      paste0(Genus, "-", row_number())
    }
  ) %>%
  ungroup() %>%
  column_to_rownames("OTU") %>%
  as.matrix()

# Add new taxonomy to phyloseq object
tax_table(ps_filt) <- tax_table(taxonomy)
rm(taxonomy)


# uses updated taxonomy to keep NA values
ps_genus = tax_glom(ps_filt, "Genus")


# ------ Save at genus level ------

# remove floccular granules
ps_sub <- subset_samples(ps_genus, size.name != "floccular")

saveRDS(ps_genus, file = "./data/phyloseq/ps_genus_full.rds")
saveRDS(ps_sub,   file = "./data/phyloseq/ps_genus_subset.rds")

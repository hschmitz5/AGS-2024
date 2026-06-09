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

ps_filt0 <- phyloseq::subset_taxa(ps,
                                 # remove unclassified sequences
                                 (Kingdom != "Unassigned") &
                                 # remove Chloroplasts and Mitochondria
                                 (Order != "Chloroplast") &
                                 (Family != "Mitochondria")
)

# Keep taxa present in at least three samples
keep_prev <- rowSums(otu_table(ps_filt0) > 0) >= 3
# Keep taxa with at least 10 total reads across all samples
keep_abund <- taxa_sums(ps_filt0) >= 10
# Apply conditions
ps_filt <- prune_taxa(keep_prev & keep_abund, ps_filt0)

# ------ Agglomerate, cutting NA values ------

# tax_glom deletes NA values by default
ps_genus_short = tax_glom(ps_filt, "Genus")

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


# ------ Save ------

saveRDS(ps_filt,        file = "./data/ps_ASV_full.rds")
saveRDS(ps_genus_short, file = "./data/ps_genus_cut_NA.rds")
saveRDS(ps_genus,       file = "./data/ps_genus_full.rds")

# remove XS granules
ps_sub <- subset_samples(ps_genus, size.name != "floccular")
  
saveRDS(ps_sub, file = "./data/ps_genus_subset.rds")

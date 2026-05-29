rm(list = ls())
library(tidyverse)
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_subset.R")

taxon <- "Ca_Contendobacter_g-18"

ps_fname    <- "./data/ps_ASV_full.rds"

# absolute counts
ps <- readRDS(ps_fname)

# Change OTU to Species_updated
taxonomy <- data.frame(tax_table(ps))
rownames(ps@otu_table) <- taxonomy$Species_updated
rownames(ps@tax_table) <- taxonomy$Species_updated
ps@phy_tree$tip.label <- taxonomy$Species_updated


# names of significant taxa
write2excel = FALSE
sig_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel, fname_excel)


ps_rel_df <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x)) %>%
  phyloseq::psmelt()

ASV_size <- ps_rel_df %>%
  dplyr::select(Sample, size.mm, size.name, OTU, Abundance) %>%
  group_by(Sample, size.mm, size.name, OTU) %>%
  summarise(
    Abundance = sum(Abundance),
    .groups = "drop"
  ) %>%
  filter(OTU %in% taxon) %>%
  group_by(size.name) %>%
  summarise(
    Avg_Abund = mean(Abundance),
    .groups = "drop"
  )

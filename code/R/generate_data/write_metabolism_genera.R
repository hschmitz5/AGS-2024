rm(list = ls())
library(writexl)
source("./code/R/01_load_ps.R")

fname_out <- "./data/midas_genera.xlsx"

# define taxa in which at least one sample has abundance > rel_ab_cutoff
high_ab_taxa <- get_rel_ASV(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

distinct_genera <- get_rel_ASV(ps) %>%
  filter(OTU %in% high_ab_taxa) %>%
  dplyr::select(Genus, OTU, Sample, Abundance) %>%  
  pivot_wider(
    names_from = Sample,
    values_from = Abundance
  ) %>%
  dplyr::select(Genus, OTU, sam_name) %>%
  filter(!is.na(Genus)) %>%
  arrange(Genus) %>%
  distinct(Genus)
 
write_xlsx(distinct_genera, path = fname_out)
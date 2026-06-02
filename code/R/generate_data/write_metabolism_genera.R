rm(list = ls())
library(writexl)
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_diff_ab.R")
source("./code/R/04_DA_agglom.R")

fname_out <- "./data/midas_genera.xlsx"

#### Rename: agglomerate names when multiple ASVs are differentially abundant or not
ancom_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold, rel_ab_cutoff, write2excel = FALSE)
DA_taxa <- ancom_taxa$high_ab 

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
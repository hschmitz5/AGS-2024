rm(list = ls())
library(writexl)
source("./code/R/01_load_ps.R")

rel_ab_cutoff <- 0.5 # percent

fname_out <- "./data/midas_genera.xlsx"

# define taxa in which at least one sample has abundance > rel_ab_cutoff
taxa_names <- get_rel(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(Genus) %>%
  pull(Genus)

high_ab_genera <- get_rel_wide(ps) %>%
  rownames_to_column(var = "Genus") %>%
  filter(Genus %in% taxa_names) %>%
  arrange(Genus) 
 
write_xlsx(high_ab_genera, path = fname_out)

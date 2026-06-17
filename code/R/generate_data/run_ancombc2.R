rm(list = ls())

library(tidyverse)
library(phyloseq)
library(ANCOMBC)

fname_out <- "./data/DA/ancombc2_genus.rds"

ps <- readRDS("./data/phyloseq/ps_genus_subset.rds") 
  
set.seed(123)
output <- ancombc2(
  data = ps,
  prv_cut = 0.1, # Removes taxa not present in 10% of samples
  fix_formula = "size.name",    
  group = "size.name",
  alpha = 0.05,
  global = TRUE, dunnet = TRUE, 
  pairwise = FALSE, trend = FALSE,
  struc_zero = FALSE
)
  
saveRDS(output, file = fname_out)


# ------ Process Data ------

taxonomy <- data.frame(tax_table(ps)) %>%
  tibble::rownames_to_column("OTU") %>%
  dplyr::select(OTU, Genus)

res_prim = output$res_dunn %>%
  rename(OTU = taxon) %>%
  left_join(taxonomy, by = "OTU")
  
df_lfc <- res_prim %>%
  # filter for values that are differentially abundant (TRUE)
  dplyr::filter(diff_size.nameM == 1 | diff_size.nameL == 1 |
                  diff_size.nameXL == 1 | diff_size.nameXXL == 1) %>%
  # set non-DA values to zero
  dplyr::mutate(M = ifelse(diff_size.nameM == 1, 
                            round(lfc_size.nameM, 2), 0),
                L = ifelse(diff_size.nameL == 1, 
                            round(lfc_size.nameL, 2), 0),
                XL = ifelse(diff_size.nameXL == 1, 
                            round(lfc_size.nameXL, 2), 0),
                XXL = ifelse(diff_size.nameXXL == 1, 
                              round(lfc_size.nameXXL, 2), 0)) %>%
  dplyr::select(Genus, M, L, XL, XXL) %>%
  tidyr::pivot_longer(!Genus, names_to = "size", values_to = "lfc") %>%
  dplyr::arrange(Genus)
  
df_robust <- res_prim %>%
  # filter for values that are differentially abundant (TRUE)
  dplyr::filter(diff_size.nameM == 1 | diff_size.nameL == 1 |
                  diff_size.nameXL == 1 | diff_size.nameXXL == 1) %>%
  # true if robust (is differentially abundant and passed sensitivity analysis)
  dplyr::mutate(M = ifelse(diff_robust_size.nameM == 1, 
                           TRUE, FALSE),
                L = ifelse(diff_robust_size.nameL == 1, 
                           TRUE, FALSE),
                XL = ifelse(diff_robust_size.nameXL == 1, 
                            TRUE, FALSE),
                XXL = ifelse(diff_robust_size.nameXXL == 1, 
                             TRUE, FALSE)) %>%
  dplyr::select(Genus, M, L, XL, XXL) %>%
  tidyr::pivot_longer(!Genus, names_to = "size", values_to = "robust") %>%
  dplyr::arrange(Genus)
  
df = df_lfc %>%
  dplyr::full_join(df_robust, by = c("Genus", "size")) %>%
  mutate(size = factor(size, levels = c("M", "L", "XL", "XXL"))) 


saveRDS(df, file = "./data/DA/DA_genus_processed.rds")

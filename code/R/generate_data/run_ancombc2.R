rm(list = ls())

library(phyloseq)
library(ANCOMBC)

fname_P <- c("./data/ps_metab_P.rds", "./data/ancombc2_metab_P.rds")
fname_PV <- c("./data/ps_metab_PV.rds", "./data/ancombc2_metab_PV.rds")

run_ancom <- function(fname_in, fname_out) {
  ps <- readRDS(fname_in) 
  
  # Keep taxa present in at least three samples
  keep_prev <- rowSums(otu_table(ps) > 0) >= 3
  # Keep taxa with at least 10 total reads across all samples
  keep_abund <- taxa_sums(ps) >= 10
  # Apply conditions
  ps_filt <- prune_taxa(keep_prev & keep_abund, ps)
  
  set.seed(123)
  output <- ancombc2(
    data = ps_filt,
    fix_formula = "size.name",    
    group = "size.name",
    struc_zero = TRUE,
    global = TRUE, dunnet = TRUE, 
    pairwise = FALSE, trend = FALSE   
  )
  
  saveRDS(output, file = fname_out)
  
  return(output)
}

output_P <- run_ancom(fname_P[1], fname_P[2])
output_PV <- run_ancom(fname_PV[1], fname_PV[2])

# ------ Process Data ------

process_DA <- function(output, value_col) {
  res_prim = output$res_dunn %>% # or res
    rename(metab = taxon)
  
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
    dplyr::select(metab, M, L, XL, XXL) %>%
    tidyr::pivot_longer(!metab, names_to = "size", values_to = "lfc") %>%
    dplyr::arrange(metab)
  
  df_robust <- res_prim %>%
    # filter for values that are differentially abundant (TRUE)
    dplyr::filter(diff_size.nameM == 1 | diff_size.nameL == 1 |
                    diff_size.nameXL == 1 | diff_size.nameXXL == 1) %>%
    # set fontface to bold for values that passed sensitivity analysis
    dplyr::mutate(M = ifelse(diff_robust_size.nameM == 1, 
                             TRUE, FALSE),
                  L = ifelse(diff_robust_size.nameL == 1, 
                             TRUE, FALSE),
                  XL = ifelse(diff_robust_size.nameXL == 1, 
                              TRUE, FALSE),
                  XXL = ifelse(diff_robust_size.nameXXL == 1, 
                               TRUE, FALSE)) %>%
    dplyr::select(metab, M, L, XL, XXL) %>%
    tidyr::pivot_longer(!metab, names_to = "size", values_to = "robust") %>%
    dplyr::arrange(metab)
  
  df = df_lfc %>%
    dplyr::full_join(df_robust, by = c("metab", "size")) %>%
    mutate(size = factor(size, levels = c("M", "L", "XL", "XXL"))) %>%
    mutate(metab_val = value_col)
  
  return(df)
}

lfc_P <- process_DA(output_P, "Positive") 
lfc_PV <- process_DA(output_PV, "Positive + Variable")
lfc_full <- bind_rows(lfc_P, lfc_PV) 

saveRDS(lfc_full, file = "./data/DA_metab_processed.rds")

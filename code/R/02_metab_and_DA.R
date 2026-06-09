library(readxl)

# ANCOM-BC2 differential abundance data
ancom_fname <- "./data/ancombc2_genus.rds"
# Metabolism input file
metab_fname <- "./data/metabolism_midas.xlsx"

# p-value used for defining significant taxa in differential abundance
p_threshold <- 0.05

# input df must contain column named Genus
get_metabolism <- function(df, metab_fname) {
  m <- read_excel(metab_fname, sheet = "input") # tibble
  
  metab <- df %>%
    dplyr::select(Genus) %>%
    distinct() %>%
    left_join(., m, by = "Genus") %>%
    column_to_rownames("Genus")
}

get_ancom_taxa <- function(fname_in, ps, p_threshold) {
  # load differential abundance data
  output <- readRDS(fname_in)
  
  taxonomy <- get_taxonomy(ps)
  
  # Subset based on significance and passing the sensitivity analysis
  all_sig_taxa <- output$res %>%
    rename(OTU = taxon) %>%
    left_join(taxonomy, join_by(OTU)) %>%
    # Combines diff_size* and passed_ss* together
    pivot_longer(
      cols = matches("q_size\\.name|passed_ss_size\\.name"),
      names_to = c(".value","size"),
      names_pattern = "(q|passed_ss)_size\\.name(.*)"
    ) %>%
    filter(q < p_threshold & passed_ss == TRUE) %>%   
    distinct(Genus) %>%
    pull(Genus)
}

# Only the difference between bias corrected data is meaningful.
# It is not a substitute for relative abundance.
get_bc_abund <- function(fname) {
  
  output <- readRDS(fname)
  
  metadata <- get_metadata(ps)
  
  bc_long <- output$bias_correct_log_table %>%
    rownames_to_column("OTU") %>%
    filter(!is.na(OTU)) %>%
    pivot_longer(-OTU, names_to = "Sample", values_to = "bc_abund") %>%
    left_join(metadata, by = "Sample") %>%
    dplyr::select(Sample, size.mm, size.name, OTU, bc_abund) 
  
  return(bc_long)
}
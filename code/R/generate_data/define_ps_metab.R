rm(list = ls())
library(readxl)
library(tidyverse)
library(phyloseq)

# For ordering
sam_name <- c("20A", "20B", "20C", "14A", "14B", "14C", "10A", "10B", "10C",
              "7A", "7B", "7C", "5A", "5B", "5C")

# load phyloseq object (absolute counts)
ps <- readRDS("./data/ps_genus_subset.rds")

taxonomy <- data.frame(tax_table(ps)) %>%
  dplyr::select(Genus) %>%
  rownames_to_column("OTU")

otu_df <- as.data.frame(as.matrix(otu_table(ps))) %>%
  rownames_to_column("OTU") %>%
  left_join(., taxonomy, by = "OTU") %>%
  dplyr::select(Genus, all_of(sam_name))

# Metabolism input file
metab_fname <- "./data/metabolism_midas.xlsx"

m <- read_excel(metab_fname, sheet = "input") %>%
  column_to_rownames("Genus")

# define taxa in each metabolism group
taxa_P <- map(m, ~ rownames(m)[which(.x == "P")])
taxa_V <- map(m, ~ rownames(m)[which(.x == "V")])
# combine taxa if the metab has entries in V
taxa_PV <- map2(taxa_P, taxa_V, ~ {
  if (length(.y) > 0) {
    union(.x, .y)   # taxa that are P or V
  } else {
    character(0)    
  }
})

metab_mat <- purrr::map_dfr(
  taxa_P,
  ~ otu_df %>%
    dplyr::filter(Genus %in% .x) %>%
    dplyr::select(-Genus) %>%
    colSums() %>%
    t() %>%
    as.data.frame(),
  .id = "metab"
) %>%
  tibble::column_to_rownames("metab")


new_otu_table <- otu_table(as.matrix(metab_mat), taxa_are_rows = TRUE)
ps_metab <- phyloseq(new_otu_table, sample_data(ps))

saveRDS(ps_metab, file = "./data/ps_metab.rds")

rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV_rarefied.rds")

# UniFrac expects a binary tree (each node to have 2 descendants)
if (ape::is.binary(phy_tree(ps)) == "FALSE") {
  phy_tree(ps) <- ape::multi2di(phy_tree(ps))
  print("resolved polytomies")
}

sizes <- c("floccular", "S", "M", "L", "XL", "XXL")

all_combos <- combn(sizes, 2, simplify = FALSE)

num_rows <- length(all_combos)

df <- data.frame(matrix(NA, nrow = num_rows, ncol = 4))
colnames(df) <- c("sz_1", "sz_2", "R2", "p_value")

for (i in seq_along(all_combos)) {
  
  sample_pair <- all_combos[[i]]
  
  ps_sub <- subset_samples(ps, size.name %in% sample_pair)
  # prune OTUs that do not exist in this pair
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  meta_sub <- data.frame(sample_data(ps_sub))

  dist_matrix <- phyloseq::distance(ps_sub, method = "unifrac", weighted = TRUE)
  
  res <- adonis2(
    dist_matrix ~ size.name,
    data = meta_sub,
    permutations = 719
  )

  df[i, ] <- tibble(
    sz_1 = sample_pair[1],
    sz_2 = sample_pair[2],
    R2 = res$R2[1],
    p_value = res$`Pr(>F)`[1]
  )
}


df_R2 <- df %>%
  dplyr::select(sz_1, sz_2, R2) %>%
  pivot_wider(
    names_from = sz_2,
    values_from = R2
  ) 

fname_out <- "./data/ADONIS_UniFrac.xlsx"
write_xlsx(df_R2, path = fname_out)

# ------ Overall ADONIS ------

metadata <- data.frame(sample_data(ps))

dist_matrix <- phyloseq::distance(ps, method = "unifrac", weighted = TRUE)

res <- adonis2(
  dist_matrix ~ size.name,
  data = metadata,
  permutations = 999
)

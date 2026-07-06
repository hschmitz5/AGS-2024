rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV.rds")

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps))  # lowest number of ASVs per sample

# ------ Overall Result ------

metadata <- data.frame(sample_data(ps))
# rows are samples, columns are OTUs
otu_matrix <- data.frame(t(otu_table(ps)))

set.seed(1)
dist_matrix <- avgdist(otu_matrix, sample = rarefy_level, iterations = 10, dmethod = "bray")

overall_res <- adonis2(
  dist_matrix ~ size.name,
  data = metadata,
  permutations = 999
  ) %>%
  rownames_to_column(var = "Data")

# Want Insignificant 
overall_bd <- anova(
  betadisper(dist_matrix, metadata$size.name)
  )

# ------ Pairwise Results ------

# Define all pairwise combinations
sizes <- c("Floccular", "S", "M", "L", "XL", "XXL")
all_combos <- combn(sizes, 2, simplify = FALSE)

# Initialize empty data frame
num_rows <- length(all_combos)
df <- data.frame(matrix(NA, nrow = num_rows, ncol = 5))
colnames(df) <- c("sz_1", "sz_2", "R2", "p_value", "bd_pval")

for (i in seq_along(all_combos)) {
  
  sample_pair <- all_combos[[i]]
  
  ps_sub <- subset_samples(ps, size.name %in% sample_pair)
  # prune OTUs that do not exist in this pair
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  meta_sub <- data.frame(sample_data(ps_sub))
  
  # rows are samples, columns are OTUs
  otu_matrix <- data.frame(t(otu_table(ps_sub)))
  
  dist_matrix <- avgdist(otu_matrix, sample = rarefy_level, iterations = 10, dmethod = "bray")
  
  res <- adonis2(
    dist_matrix ~ size.name,
    data = meta_sub,
    permutations = 719
  )
  
  # Want Insignificant 
  bd_res <- anova(
    betadisper(dist_matrix, meta_sub$size.name)
    )

  df[i, ] <- tibble(
    sz_1 = sample_pair[1],
    sz_2 = sample_pair[2],
    R2 = res$R2[1],
    p_value = res$`Pr(>F)`[1],
    bd_pval = bd_res$'Pr(>F)'[1]
  )
}


df_p <- df %>%
  dplyr::select(sz_1, sz_2, p_value) %>%
  pivot_wider(
    names_from = sz_2,
    values_from = p_value
  ) 

df_R2 <- df %>%
  dplyr::select(sz_1, sz_2, R2) %>%
  pivot_wider(
    names_from = sz_2,
    values_from = R2
  ) 

df_bd <- df %>%
  dplyr::select(sz_1, sz_2, bd_pval) %>%
  pivot_wider(
    names_from = sz_2,
    values_from = bd_pval
  ) 

saveRDS(df_bd, file = "./data/bray_betadisper_result.rds")

bd_long <- df_bd |>
  pivot_longer(!sz_1, names_to = "sz_2", values_to = "bd") 

bd_long$sz_1 <- factor(bd_long$sz_1, levels = rev(sizes))
bd_long$sz_2 <- factor(bd_long$sz_2, levels = sizes)

# ------ Write to Excel ------

fname_out <- "./data/ADONIS_Bray.xlsx"
write_xlsx(
  list(
    p_values = df_p,
    bd_pvalues = df_bd,
    R2 = df_R2,
    overall = overall_res,
    bd = overall_bd
  ),
  path = fname_out
  )

# ------ Plot ------

p <- ggplot(data = bd_long, aes(x = sz_2, y = sz_1, fill = bd)) +
  geom_tile() + 
  geom_text(aes(label = round(bd,2))) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    na.value = "white"
  ) +
  labs(
    title = "Dispersion Homogeneity p-values",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

fname <- "./figures/betadisper_Bray.png"
ggsave(fname, plot = p, width = 5.38, height = 3, dpi = 300)

# ------ Look at distance matrix ------

# For ordering
sam_name <- c("20A", "20B", "20C", "14A", "14B", "14C", "10A", "10B", "10C",
              "7A", "7B", "7C", "5A", "5B", "5C")

dist_mat <- as.matrix(dist_matrix)
dist_mat <- dist_mat[sam_name, sam_name] # order data
dist_tbl <- as_tibble(dist_mat, rownames = "sample") 

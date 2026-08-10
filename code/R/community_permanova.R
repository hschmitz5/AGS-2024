rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)

# ------ Overall Result ------

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV.rds")

metadata <- data.frame(sample_data(ps))

# rows are samples, columns are OTUs
otu_matrix_full <- t(
  as.data.frame(otu_table(ps))
)

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps))  # lowest number of ASVs per sample

set.seed(1)
dist_matrix_full <- avgdist(otu_matrix_full, sample = rarefy_level, iterations = 10, dmethod = "bray")

# PERMANOVA
overall_res <- adonis2(
  dist_matrix_full ~ size.name,
  data = metadata,
  permutations = 999
  ) %>%
  rownames_to_column(var = "Data")

# Multivariate homogeneity of groups dispersions
overall_bd <- anova(
  betadisper(dist_matrix_full, metadata$size.name)
  )

# ------ Pairwise Results ------

# Define all pairwise combinations
sizes <- c("Floccular", "S", "M", "L", "XL", "XXL")
all_combos <- combn(sizes, 2, simplify = FALSE)

pairwise_adonis <- function(sample_pair) {
  
  keep <- sample_names(ps)[sample_data(ps)$size.name %in% sample_pair]
  
  ps_sub <- prune_samples(keep, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  
  meta <- data.frame(sample_data(ps_sub))
  otu <- t(as.data.frame(otu_table(ps_sub)))
  
  dist <- avgdist(
    otu,
    sample = rarefy_level,
    iterations = 10,
    dmethod = "bray"
  )
  
  ad <- adonis2(
    dist ~ size.name,
    data = meta,
    permutations = 719
  )
  
  bd <- anova(
    betadisper(dist, meta$size.name)
  )
  
  tibble(
    sz_1    = sample_pair[1],
    sz_2    = sample_pair[2],
    R2      = ad$R2[1],
    p_value = ad$`Pr(>F)`[1],
    bd_pval = bd$`Pr(>F)`[1]
  )
}

df_long <- map_dfr(all_combos, pairwise_adonis) %>%
  mutate(
    sz_1 = factor(sz_1, levels = rev(sizes)),
    sz_2 = factor(sz_2, levels = sizes),
  )

# Create list that combines pairwise statistics
stat_var <- c("p_value", "R2", "bd_pval")

df_wide <- map(
  set_names(stat_var), \(var) {
  df_long %>%
    dplyr::select(sz_1, sz_2, all_of(var)) %>%
    pivot_wider(
      names_from = sz_2,
      values_from = all_of(var)
      )
    })

df_wide$p_value
df_wide$R2
df_wide$bd_pval

# ------ Plot ------

bd_long <- df_long %>%
  select(sz_1, sz_2, bd = bd_pval) 

p <- ggplot(data = bd_long, aes(x = sz_2, y = sz_1, fill = bd)) +
  geom_tile() + 
  geom_text(aes(label = round(bd,2))) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    na.value = "white"
  ) +
  labs(
    title = expression("Dispersion Homogeneity " * italic(p) * "-values"), #"Dispersion Homogeneity p-values",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

fname <- "./figures/Figure_S.2.tif"
ggsave(fname, plot = p, width = 5.38, height = 3, dpi = 300)


# ------ Write to Excel ------

fname_out <- "./data/ADONIS_Bray.xlsx"
write_xlsx(
  list(
    overall = overall_res,
    bd = overall_bd,
    p_values = df_wide$p_value,
    R2 = df_wide$R2,
    bd_pvalues = df_wide$bd_pval
  ),
  path = fname_out
)

# ------ Look at distance matrix ------

# For ordering
sam_name <- c("20A", "20B", "20C", "14A", "14B", "14C", "10A", "10B", "10C",
              "7A", "7B", "7C", "5A", "5B", "5C")

dist_mat <- as.matrix(dist_matrix_full)
dist_mat <- dist_mat[sam_name, sam_name] # order data
dist_tbl <- as_tibble(dist_mat, rownames = "sample") 

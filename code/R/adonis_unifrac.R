rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV_rarefied.rds")

# ------ Overall Result ------

metadata <- data.frame(sample_data(ps))

dist_matrix <- phyloseq::distance(ps, method = "wunifrac") 

overall_res <- adonis2(
  dist_matrix ~ size.name,
  data = metadata,
  permutations = 999
)

bd <- betadisper(dist_matrix, metadata$size.name)
overall_bd <- anova(bd)

# ------ Pairwise Results ------

# Define all pairwise combinations
sizes <- c("Floccular", "S", "M", "L", "XL", "XXL")
all_combos <- combn(sizes, 2, simplify = FALSE)

# Initialize empty data frame
num_rows <- length(all_combos)
df <- data.frame(matrix(NA, nrow = num_rows, ncol = 4))
colnames(df) <- c("sz_1", "sz_2", "R2", "p_value")

for (i in seq_along(all_combos)) {
  
  sample_pair <- all_combos[[i]]
  
  ps_sub <- subset_samples(ps, size.name %in% sample_pair)
  # prune OTUs that do not exist in this pair
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  meta_sub <- data.frame(sample_data(ps_sub))

  dist_matrix <- phyloseq::distance(ps_sub, method = "wunifrac")
  
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

# ------ Write to Excel ------

fname_out <- "./data/ADONIS_UniFrac.xlsx"
write_xlsx(
  list(
    p_values = df_p,
    R2 = df_R2,
    overall = overall_res,
    bd = overall_bd
  ),
  path = fname_out
)

# ------ Plot ------

df_long <- df_p |>
  pivot_longer(!sz_1, names_to = "sz_2", values_to = "p_val") %>%
  mutate(sz_1 = recode(sz_1, "floccular" = "Floccular"))

df_long$sz_1 <- factor(df_long$sz_1, levels = rev(sizes))
df_long$sz_2 <- factor(df_long$sz_2, levels = sizes)
  
p <- ggplot(data = df_long, aes(x = sz_2, y = sz_1, fill = p_val)) +
  geom_tile() + 
  geom_text(aes(label = p_val)) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    na.value = "white"
  ) +
  labs(
    title = "Permutational MANOVA p-values",
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

fname <- "./figures/ADONIS_wUniFrac_p.png"
ggsave(fname, plot = p, width = 6.5, height = 3, dpi = 300)

rm(list = ls())
library(vegan)

# load phyloseq object for all sample sizes
ps_full <- readRDS("./data/ps_ASV_full.rds") 

metadata <- data.frame(sample_data(ps_full)) %>%
  tibble::rownames_to_column("Sample") %>%
  arrange(size.name)

ps.dist <- distance(ps_full, method = "wunifrac" ) # weighted
invisible(ps.dist)
#show(ps.dist)

ps.ano <- anosim(ps.dist, metadata$size.name, permutations = 9999)

# Get summary and write to a text file
txt_fname <- "./figures/ANOSIM.txt"
anosim_summary <- capture.output(summary(ps.ano))
writeLines(anosim_summary, con = txt_fname)

# ------- Plot --------

# Define a base font size in points
text_size <- 10

# Extract ranks and groups into a data frame
ano_df <- data.frame(
  rank = ps.ano$dis.rank,           # ranked dissimilarities
  group = ps.ano$class.vec          # group labels
)

anosim_plot <- ggplot(ano_df, aes(x = group, y = rank)) +
  # Boxplot for Between points
  geom_boxplot(data = subset(ano_df, group == "Between"),
               aes(x = group, y = rank),
               alpha = 0.3, outlier.shape = NA) +
  # Scatter for non-Between points
  geom_point(data = subset(ano_df, group != "Between"),
              width = 0.1, height = 0) +
  labs(
    title = paste0("R = ", round(ps.ano$statistic, 3), ", P = ", ps.ano$signif),
    x = "Size",
    y = "Dissimilarity Ranks"
  ) +
  theme_minimal(base_size = text_size)

ano_fname <- "./figures/anosim_plot.png"
ggsave(ano_fname, plot = anosim_plot, width = 4, height = 3, dpi = 300)

rm(list = ls())
source("./code/R/01_load_ps.R")

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV_rarefied.rds") 

# ps: all sample groups
ps.ord <- ordinate(ps, "PCoA", "wunifrac")

# symbol
# 16 = filled circle, 17 = triangle, 15 = square, 18 = diamond, etc.
shapes <- c(16, 17, 15, 18, 3, 7)

# colors
cols <- c("gray", met.brewer(size_pal, n_sizes))

p <- plot_ordination(ps, ps.ord, shape = "size.name", color="size.name") +
  geom_polygon(alpha = 0.5, aes(fill = size.name)) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(shape = "Size", color = "Size", fill = "Size") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12)
  ) 

fname_ord <- "./figures/ordination-PCoA-ASV.png"
ggsave(fname_ord, plot = p, width = 6.5, height = 3, dpi = 600)


# For correlation
metadata <- get_metadata(ps) 

pcoa <- data.frame(axis1 = ps.ord$vectors[, "Axis.1"]) %>%
  rownames_to_column(var = "Sample") %>%
  left_join(., metadata, by = "Sample") %>%
  arrange(size.mm)

res <- cor.test(pcoa$axis1, pcoa$size.midpoint, method = "spearman")

print(res$p.value)

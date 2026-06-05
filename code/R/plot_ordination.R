rm(list = ls())
source("./code/R/01_load_ps.R")

fname_ord <- "./figures/ordination-PCoA.png"
fname_ord2 <- "./figures/ordination-PCoA-mu.png"

# load phyloseq object for all sample sizes
ps_full <- readRDS("./data/ps_ASV_full.rds")

# ps_full: all sample groups
ps.ord <- ordinate(ps_full, "PCoA", "wunifrac")

# colors
cols <- c("gray", met.brewer(size_pal, n_sizes))

# symbol
# 16 = filled circle, 17 = triangle, 15 = square, 18 = diamond, etc.
shapes <- c(16, 17, 15, 18, 3, 7)

p <- plot_ordination(ps_full, ps.ord, type="samples", color="size.name", shape = "size.name") +
  # stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=size.name)) + # need at least 4 points
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  labs(color = "Size", shape = "Size") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 12)
  ) 

ordination_plot <- p

ggsave(fname_ord, plot = ordination_plot, width = 6.5, height = 3, dpi = 600)


# For correlation
metadata <- get_metadata(ps_full) 

pcoa <- data.frame(axis1 = ps.ord$vectors[, "Axis.1"]) %>%
  rownames_to_column(var = "Sample") %>%
  left_join(., metadata, by = "Sample") %>%
  arrange(size.mm)

res <- cor.test(pcoa$axis1, pcoa$size.midpoint, method = "spearman")

print(res$p.value)

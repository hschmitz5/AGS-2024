rm(list = ls())
source("./code/R/01_load_ps.R")

dist_method <- "bray" # "wunifrac"

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV_rarefied.rds") 

# UniFrac expects a binary tree (each node to have 2 descendants)
# if (ape::is.binary(phy_tree(ps)) == "FALSE") {
#   phy_tree(ps) <- ape::multi2di(phy_tree(ps))
#   print("resolved polytomies")
# }

# ps: all sample groups
ps.ord <- ordinate(ps, "PCoA", "bray") 

# ------ Correlation ------

# For correlation
metadata <- get_metadata(ps) 

pcoa <- data.frame(axis1 = ps.ord$vectors[, "Axis.1"]) %>%
  rownames_to_column(var = "Sample") %>%
  left_join(., metadata, by = "Sample") %>%
  group_by(size.name, size.midpoint) %>%
  # correlate mean of Axis 1
  summarize(
    axis1_mean = mean(axis1),
    .groups = "drop"
  )

res <- cor.test(pcoa$axis1_mean, pcoa$size.midpoint, method = "spearman")

print(res$p.value)

# ------ Plot ------

# symbol
# 16 = filled circle, 17 = triangle, 15 = square, 18 = diamond, etc.
shapes <- c(16, 17, 15, 18, 3, 7)

# colors
cols <- c("gray", met.brewer(size_pal, n_sizes))

p <- plot_ordination(ps, ps.ord, color="size.name", shape = "size.name") +
  geom_polygon(alpha = 0.5, aes(fill = size.name)) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = cols) +
  labs(shape = "Size", color = "Size", fill = "Size") +
  theme_classic(base_size = 12) 

fname_ord <- paste0("./figures/ordination_",dist_method,".png")
ggsave(fname_ord, plot = p, width = 6.5, height = 3, dpi = 300)

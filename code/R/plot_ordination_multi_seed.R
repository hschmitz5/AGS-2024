rm(list = ls())
source("./code/R/01_load_ps.R")

# load phyloseq object for all sample sizes
ps_asv <- readRDS("./data/phyloseq/ps_ASV.rds") 

# ------ Rarefy x 10 ------

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps_asv))  # lowest number of ASVs per sample

seeds <- set_names(1:10, paste0("seed_", 1:10))

# apply rarefaction ten times
ps_rar_list <- map(seeds, ~ rarefy_even_depth(ps_asv, sample.size = rarefy_level,
                                             rngseed = .x, replace = FALSE, trimOTUs = FALSE, verbose = TRUE))


# ------ Compute distance per rarefaction -----

for (seed_num in 1:10) {
  fname_ord <- paste0("./figures/ordination by seed/ordination-PCoA-ASV-seed",seed_num,".png")
  
  ps_seed <- ps_rar_list[[paste0("seed_", seed_num)]]
  
  # ps_full: all sample groups
  ps.ord <- ordinate(ps_seed, "PCoA", "wunifrac")
  
  # symbol
  # 16 = filled circle, 17 = triangle, 115 = square, 18 = diamond, etc.
  shapes <- c(16, 17, 15, 18, 3, 7)
  
  # colors
  cols <- c("gray", met.brewer(size_pal, n_sizes))
  
  p <- plot_ordination(ps_seed, ps.ord, shape = "size.name", color="size.name") +
    geom_polygon(alpha = 0.5, aes(fill = size.name)) +
    scale_shape_manual(values = shapes) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(shape = "Size", color = "Size", fill = "Size") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 12)
    ) 
  
  ordination_plot <- p 
  
  ggsave(fname_ord, plot = ordination_plot, width = 6.5, height = 3, dpi = 600)
}

# ps_seed <- ps_rar_list[["seed_1"]]
# 
# # For correlation
# metadata <- get_metadata(ps_seed) 
# 
# pcoa <- data.frame(axis1 = ps.ord$vectors[, "Axis.1"]) %>%
#   rownames_to_column(var = "Sample") %>%
#   left_join(., metadata, by = "Sample") %>%
#   arrange(size.mm)
# 
# res <- cor.test(pcoa$axis1, pcoa$size.midpoint, method = "spearman")
# 
# print(res$p.value)

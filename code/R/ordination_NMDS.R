rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)
source("./code/R/01_load_ps.R")

# load phyloseq object for all sample sizes
ps <- readRDS("./data/phyloseq/ps_ASV.rds")

metadata <- data.frame(sample_data(ps)) %>%
  rownames_to_column(var = "SampleID")

# rows are samples, columns are OTUs
otu_matrix <- data.frame(t(otu_table(ps)))

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps))  # lowest number of ASVs per sample

set.seed(1)
dist_matrix <- avgdist(otu_matrix, sample = rarefy_level, iterations = 10, dmethod = "bray")

set.seed(2)
nmds <- metaMDS(dist_matrix) # list

nmds_df <- scores(nmds) %>%
  as_tibble(rownames = "SampleID") %>%
  left_join(., metadata, by = "SampleID")

stress_text <- paste0("2D stress: ",round(nmds$stress,2))

# ------ Plot ------

# symbol
# 16 = filled circle, 17 = triangle, 15 = square, 18 = diamond, etc.
shapes <- c(16, 17, 15, 18, 3, 7)

# colors
cols <- c("gray", met.brewer(size_pal, n_sizes))

p <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = size.name, shape = size.name)) +
  geom_point() +
  geom_polygon(alpha = 0.5, aes(fill = size.name)) +
  # annotate("text", x = Inf, y = -Inf, label = stress_text,
  #          hjust = 1.1, vjust = -0.5) +
  scale_color_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Axis 1", 
    y = "Axis 2",
    shape = "Size", 
    color = "Size", 
    fill = "Size") +
  theme_classic(base_size = 12) 

fname <- "./figures/ordination_NMDS.png"
ggsave(fname, plot = p, width = 6.5, height = 3, dpi = 300)

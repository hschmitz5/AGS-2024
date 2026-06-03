suppressPackageStartupMessages({
  # ---- Load Packages ----
  library(tidyverse)
  library(phyloseq)
  # formatting figures
  library(patchwork)
  # colors
  library(RColorBrewer)
  library(MetBrewer) 
})
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache = TRUE)

# Color Palettes (MetBrewer)
size_pal <- "Greek"     # reds
taxa_pal <- "Hiroshige" # orange, blue

# Relative Abundance Cutoff (%) 
rel_ab_cutoff <- 0.5

# p-value used for defining significant taxa in differential abundance
p_threshold   <- 0.05

# load phyloseq object (absolute counts)
ps_fname <- "./data/ps_ASV_subset.rds"
ps <- readRDS(ps_fname)

# ------ Define size groups ----------

# define dimensions of sample grouping
n_replicates <- 3
n_sizes      <- length(levels(ps@sam_data$size.mm))

# define sample names
size <- data.frame(
  ranges = levels(ps@sam_data$size.mm),
  name = levels(ps@sam_data$size.name)
)

# For ordering
sam_name <- c("20A", "20B", "20C", "14A", "14B", "14C", "10A", "10B", "10C",
              "7A", "7B", "7C", "5A", "5B", "5C")

# ------ Change OTU to Species_updated --------

taxonomy <- data.frame(tax_table(ps))
rownames(ps@otu_table) <- taxonomy$Species_updated
rownames(ps@tax_table) <- taxonomy$Species_updated
ps@phy_tree$tip.label <- taxonomy$Species_updated
rm(taxonomy)

# ------ Functions for processing phyloseq object ----

get_metadata <- function(ps){
  metadata <- data.frame(sample_data(ps)) %>%
    tibble::rownames_to_column("Sample") %>%
    arrange(size.name)
}

get_taxonomy <- function(ps){
  taxonomy <- data.frame(tax_table(ps)) %>%
    tibble::rownames_to_column("OTU")
}

get_rel_ASV <- function(ps) {
  # define relative abundance
  ps_rel_ASV <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x))
  
  phyloseq::psmelt(ps_rel_ASV)
}

get_diversity <- function(ps) {
  library(picante) 
  # metadata
  metadata <- get_metadata(ps)
  # community data
  comm <- as.data.frame(as.matrix(ps@otu_table)) %>%
    t() 
  # phylogenetic tree
  phy <- ps@phy_tree
  
  # PD: phylogenetic diversity (Faith's)
  # SR: species richness
  pd_results <- pd(comm, phy) %>%
    dplyr::select(PD) %>%
    rownames_to_column(var = "Sample") %>%
    mutate(Sample = factor(Sample, levels = sam_name)) %>%
    arrange(Sample) %>%
    column_to_rownames(var = "Sample")
}

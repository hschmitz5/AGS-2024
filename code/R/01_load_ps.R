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


# load phyloseq object (absolute counts)
ps_fname <- "./data/ps_genus_subset.rds"
ps <- readRDS(ps_fname)

metadata <- data.frame(sample_data(ps)) %>%
  rownames_to_column(var = "Sample") 


# Color Palettes (MetBrewer)
size_pal <- "Greek"     # reds
taxa_pal <- "Hiroshige" # orange, blue


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


# ------ Functions for processing phyloseq object ----

get_metadata <- function(ps){
  # cannot use as.data.frame because class will be phyloseq 
  metadata <- data.frame(sample_data(ps)) %>%
    tibble::rownames_to_column("Sample") %>%
    arrange(size.name)
}

get_taxonomy <- function(ps){
  taxonomy <- data.frame(tax_table(ps)) %>%
    tibble::rownames_to_column("OTU")
}

get_rel <- function(ps) {
  # define relative abundance
  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x))
  # combines taxonomy and abundance
  phyloseq::psmelt(ps_rel)
}

get_rel_wide <- function(ps) {
  rel_wide <- get_rel(ps) %>%
    dplyr::select(Genus, Sample, Abundance) %>%  
    pivot_wider(
      names_from = Sample,
      values_from = Abundance
    ) %>%
    column_to_rownames(var = "Genus") %>%
    dplyr::select(all_of(sam_name)) 
}

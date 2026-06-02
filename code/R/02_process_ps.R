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
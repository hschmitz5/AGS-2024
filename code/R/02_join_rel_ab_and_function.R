# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

join_rel_ab_and_function <- function(ps) {
  
  # Define relative abundance
  rel_size <- get_rel(ps) 
  
  # Load metabolism data
  # Input must contain Genus
  m <- get_metabolism(rel_size)
  
  # define taxa in each metabolism group
  taxa_P <- map(m, ~ rownames(m)[which(.x == "P")])
  taxa_V <- map(m, ~ rownames(m)[which(.x == "V")])
  # combine taxa if the metab has entries in V
  taxa_PV <- map2(taxa_P, taxa_V, ~ {
    if (length(.y) > 0) {
      union(.x, .y)   # taxa that are P or V
    } else {
      character(0)    
    }
  })
  
  # function sums up values when taxa = P or V
  sum_metab <- function(taxa_list) {
    map_dfr(names(taxa_list), function(nm) {
      rel_size %>%
        filter(Genus %in% taxa_list[[nm]]) %>%
        group_by(size.name, Sample) %>%
        summarize(
          sum_abund = sum(Abundance, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(metab = nm)
    }) 
  }
  
  df_P  <- sum_metab(taxa_P)  
  # P + V: only sum metab found in V
  df_PV <- sum_metab(taxa_PV) 
  
  # Summarize for plotting
  summarize_metab <- function(df, value_col) {
    df %>%
      group_by(metab, size.name) %>%
      summarize(
        mean_sum = mean(sum_abund, na.rm = TRUE),
        sd_sum = sd(sum_abund, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(metab_val = value_col)
  }
  
  P_summary <- summarize_metab(df_P, "Positive")
  PV_summary <- summarize_metab(df_PV, "Positive + Variable")
  
  # joins data sets
  full_summary <- bind_rows(P_summary, PV_summary) %>%
    # define panel grouping
    mutate(
      size.name = factor(size.name, levels = size$name)
    ) 

}

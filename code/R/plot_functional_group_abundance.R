# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

rm(list = ls())
library(patchwork)
source("./code/R/01_load_ps.R")
source("./code/R/02_join_rel_ab_and_function.R")

write2excel <- 0

metab_order <- c("GAO", "PAO", "AOB", "NOB", "Filamentous")
n_rows <- length(metab_order)

DA_df <- readRDS("./data/DA/DA_metab_processed.rds") %>%
  mutate(metab = factor(metab, levels = metab_order))
    
rel_ab_df <- join_rel_ab_and_function(ps) %>%
  filter(metab %in% metab_order) %>%
  mutate(
    metab = factor(metab, levels = metab_order)
    )

# ------------ Plot ------------------

min_y = floor(min(DA_df$lfc))
max_y = ceiling(max(DA_df$lfc))

p1 <- ggplot(rel_ab_df, aes(x = size.name, y = mean_sum, fill = metab_val)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_sum - sd_sum, ymax = mean_sum + sd_sum),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~metab, scales = "free_y", ncol = 1) +
  labs(
    title = "Relative Abundance",
    y = "Percent of Sample",
    x = "Size"
  ) 

p2 <- ggplot(DA_df, aes(x = size, y = lfc, fill = metab_val)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "darkgray") +  # bold y = 0
  facet_wrap(~metab, scales = "fixed", ncol = 1) +
  ylim(min_y, max_y) +
  labs(
    title = "Differential Abundance",
    y = "Log Fold-Change (Relative to S)",
    x = "Size"
    ) 

p <- p1 + p2 +
  # Combine legends
  plot_layout(guides = "collect") & 
  scale_fill_manual(
    name = "Functional Group",
    values = c("Positive" = "steelblue",
               "Positive + Variable" = "lightgray")
  ) &
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5), # center title
    legend.position = "bottom",
    strip.background = element_rect(
      colour = NA # facet label outline
      )
    ) 
  
# Save plot
fname <- "./figures/functional_group_abundance.png"
ggsave(fname, plot = p, width = 6.5, height = 7.5, dpi = 300)



# ------ Write Data to Excel

if (write2excel == 1) {
  ### Do not exclude data
  # define relative abundance
  rel_wide <- get_rel_wide(ps) %>%
    rownames_to_column(var = "Genus")
  
  new_m <- get_metabolism(rel_wide) %>%
    # true if any metabolic groups in row are defined
    mutate(tf = as.integer(if_any(everything(), ~ !is.na(.x)))) %>%
    rownames_to_column(var = "Genus")
  
  full_df <- left_join(rel_wide, new_m, by = "Genus") %>%
    filter(tf == 1) %>%
    dplyr::select(-tf) %>%
    relocate(where(is.numeric), .after = where(is.character)) %>%
    arrange(Genus)
  
  library(writexl)
  write_xlsx(full_df, path = "./data/functional_rel_ab.xlsx")
}
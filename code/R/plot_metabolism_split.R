rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_metab_and_DA.R")

# Metabolism input file
metab_fname <- "./data/metabolism_midas.xlsx"

# define taxa in which at least one sample has abundance > rel_ab_cutoff
high_ab_taxa <- get_rel_ASV(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

# Define relative abundance
ASV_size <- get_rel_ASV(ps) %>%
  filter(OTU %in% high_ab_taxa) %>%
  group_by(Genus, OTU, size.name) %>%
  summarize(
    mean_ab = mean(Abundance),
    std_ab  = sd(Abundance),
    .groups = "drop"
  ) 

# Load metabolism data
# Input must contain Genus
m <- get_metabolism(ASV_size, metab_fname)

# define taxa in each metabolism group
taxa_P <- map(m, ~ rownames(m)[which(.x == "P")])
taxa_V <- map(m, ~ rownames(m)[which(.x == "V")])

# function sums up values when taxa = P or V
summarize_metab <- function(taxa_list, value_col) {
  map_dfr(names(taxa_list), function(nm) {
    ASV_size %>%
      filter(OTU %in% taxa_list[[nm]]) %>%
      group_by(size.name) %>%
      summarize(
        sum_mean = sum(mean_ab, na.rm = TRUE),
        sum_sd   = sum(std_ab, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(metab_val = value_col, metab = nm) 
  }) 
}

df_P <- summarize_metab(taxa_P, "P")
df_V <- summarize_metab(taxa_V, "V")
# joins data sets
df <- bind_rows(df_P, df_V) %>%
  # define panel grouping
  mutate(
    panel = case_when(
      metab %in% c("GAO", "Nitrite reduction") ~ "Nitrite reduction & GAO",
      metab %in% c("PAO", "Filamentous") ~ "Filamentous & PAO",
      metab %in% c("AOB", "NOB") ~ "AOB & NOB"
    ),
    metab_val = recode(metab_val, "P" = "Positive", "V" = "Variable")
  ) %>%
  # reorder
  dplyr::select(panel, metab, metab_val, size.name, sum_mean, sum_sd)

# Convert to factor
df$size.name <- factor(df$size.name, levels = size$name)
df$metab <- factor(
  df$metab, levels = c("GAO", "PAO", "Nitrite reduction", "Filamentous", "AOB", "NOB")
)
df$panel <- factor(
  df$panel, levels = c("Nitrite reduction & GAO", "Filamentous & PAO", "AOB & NOB")
)

# ------------ Plot ------------------

p <- ggplot(data = df, 
            aes(x = size.name, y = sum_mean, linetype = metab_val, group = metab_val)) +
  geom_point(color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_errorbar(
    aes(ymin = sum_mean - sum_sd, ymax = sum_mean + sum_sd),
    color = "steelblue",
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  facet_wrap(~metab, nrow = 3, scales = "free") +
  labs(
    x = "Size",
    y = "Relative Abundance [%]"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.spacing.x = unit(0.5, "in"),
    legend.title = element_blank()
  )

# Save plot
fname <- "./figures/metabolism.png"
ggsave(fname, plot = p, width = 6.5, height = 5, dpi = 300)

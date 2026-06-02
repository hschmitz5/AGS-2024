rm(list = ls())
source("./code/R/00_setup.R")
source("./code/R/01_load_data.R")
source("./code/R/02_process_ps.R")
source("./code/R/03_diff_ab.R")

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
taxa <- map(m, ~ rownames(m)[!is.na(.x)])

df <- map_dfr(names(taxa), function(nm) {
  ASV_size %>%
    filter(OTU %in% taxa[[nm]]) %>%
    group_by(size.name) %>%
    summarize(
      sum_mean = sum(mean_ab),
      sum_sd   = sum(std_ab),
      .groups = "drop"
    ) %>%
    mutate(metab = nm)
  }) %>%
  # remove AOB (only one ASV which is same as NOB)
  filter(!metab %in% "AOB") %>%
  mutate(
    metab = if_else(metab == "NOB", "AOB = NOB", metab),
    panel = case_when(
      metab %in% c("GAO", "Nitrite reduction") ~ "Nitrite reduction & GAO",
      metab %in% c("PAO", "Filamentous") ~ "Filamentous & PAO",
      metab %in% c("AOB = NOB") ~ "AOB & NOB"
    )
  )

df$metab <- factor(
  df$metab, levels = c("Nitrite reduction", "GAO", "Filamentous", "PAO", "AOB = NOB")
  )

df$panel <- factor(
  df$panel, levels = c("Nitrite reduction & GAO", "Filamentous & PAO", "AOB & NOB")
  )

# ------------ Plot ------------------

p <- ggplot(data = df,
            aes(x = size.name, y = sum_mean, color = metab, group = metab)) +
  geom_point() +
  geom_line() +
  geom_errorbar(
    aes(ymin = sum_mean - sum_sd, ymax = sum_mean + sum_sd),
    width = 0.2,
    position = position_dodge(width = 0.2)
  ) +
  facet_wrap(~panel, nrow = 2, scales = "free") +
  labs(
    x = "Size",
    y = "Relative Abundance [%]"
  ) +
  scale_color_manual(
    values = c(
      "Nitrite reduction" = "black", "GAO" = "rosybrown",
      "Filamentous" = "gray", "PAO" = "darkorchid4", 
      "AOB = NOB" = "lightsteelblue"
    )
  ) +
  guides(color = guide_legend(nrow = 2, byrow = FALSE)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.spacing.x = unit(0.5, "in"),
    legend.title = element_blank()
  ) 

# Save plot
fname <- "./figures/metabolism.png"
ggsave(fname, plot = p, width = 6.5, height = 5, dpi = 300)

rm(list = ls())
library(tidyverse)
library(readxl)
library(cowplot)
library(MetBrewer)
source("./code/R/01_load_ps.R")
rm(ps)

fname_in <- "./data/Rheometry_Nov_2024.xlsx"

size_meta <- size_meta %>%
  rename(size = name)

modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  left_join(size_meta, by = "size") %>%
  select(midpoint, size, freq_rad, G_1, G_2, G_3, G2_1, G2_2, G2_3) %>%
  filter(size != "Floccular") %>%
  pivot_longer(
    cols = c(G_1:G_3, G2_1:G2_3),
    names_to = c(".value", "replicate"),
    names_pattern = "(G2?)_(\\d)"
  ) 

# ------ Summarize mean across replicates for plotting ------

mod_summary_long <- modulus %>%
  group_by(midpoint, size, freq_rad) %>%
  summarize(
    G_avg = mean(G),
    G_sd = sd(G),
    G2_avg = mean(G2),
    G2_sd = sd(G2),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(G_avg, G_sd, G2_avg, G2_sd),
    names_to = c("measure", ".value"),
    names_pattern = "(G2?|G2?)_(avg|sd)"
  ) %>%
  mutate(
    # convert units to kPa (originally in Pa)
    avg = avg/1000, 
    sd = sd/1000,
    # change display names and order
    size = factor(size, levels = size_meta$size),
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus (G')","G2"='Loss Modulus (G")')
  )

#### Plot

p1 <- ggplot(mod_summary_long, aes(x = freq_rad, y = avg, color = size)) +
  geom_point() +
  geom_line(aes(group = size)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_color_manual(
    name = "Size", 
    values = met.brewer(size_pal, n_sizes)
  ) +
  labs(
    x = "Frequency (rad/s)",
    y = "Modulus (kPa)",
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      colour = NA # facet label outline
    ),
    legend.position = "none")

fname_out <- "./figures/moduli_small.tif"
ggsave(fname_out, plot = p1, width = 4, height = 2, dpi = 300)

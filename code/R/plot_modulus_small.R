rm(list = ls())
library(tidyverse)
library(readxl)
library(cowplot)
library(MetBrewer)
source("./code/R/01_load_ps.R")
rm(ps)

fname_in <- "./data/Rheometry_Nov_2024.xlsx"

modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  filter(size != "Floccular") %>%
  select(size, freq_rad, G_avg, G_sd, G2_avg, G2_sd) %>%
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
    size = factor(size, levels = size_meta$name),
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus (G')","G2"='Loss Modulus (G")')
  )

modulus_subset <- modulus %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

#### Plot

p1 <- ggplot(modulus, aes(x = freq_rad, y = avg, color = size)) +
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

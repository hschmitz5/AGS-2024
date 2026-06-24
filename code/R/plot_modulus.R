rm(list = ls())
library(readxl)
library(tidyverse)
library(patchwork)
library(MetBrewer)
source("./code/R/01_load_ps.R")

# define sample names
sz <- data.frame(
  name = c("S", "M", "L", "XL", "XXL")
)

fname_in <- "./data/Rheometry_Nov_2024.xlsx"
modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  filter(size != "XS") %>%
  select(size, freq_rad, G_avg, G_sd, G2_avg, G2_sd) %>%
  pivot_longer(
    cols = c(G_avg, G_sd, G2_avg, G2_sd),
    names_to = c("measure", ".value"),
    names_pattern = "(G2?|G2?)_(avg|sd)"
  ) %>%
  mutate(
    size = factor(size, levels = sz$name),
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus (G')","G2"='Loss Modulus (G")')
    )

modulus_subset <- modulus %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

#### Plot

p <- ggplot(modulus, aes(x = freq_rad, y = avg, color = size)) +
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
    y = "Modulus (Pa)",
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right")

fname_out <- "./figures/moduli.png"
ggsave(fname_out, plot = p, width = 6.5, height = 2.25, dpi = 300)


p_sub <- ggplot(modulus_subset, aes(x = size, y = avg, fill = measure)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  labs(
    title = "Frequency = 0.1 rad/s",
    x = "Size",
    y = "Modulus (Pa)"
  ) +
  scale_fill_manual(
    values = c("#765DA0", "lightgray")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

fname_out <- "./figures/moduli_subset.png"
ggsave(fname_out, plot = p_sub, width = 5, height = 2.25, dpi = 300)

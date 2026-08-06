rm(list = ls())
library(tidyverse)
library(readxl)
library(cowplot)
library(MetBrewer)
source("./code/R/01_load_ps.R")
rm(ps)

raw_df <- read_excel("./data/Rheometry_Nov_2024.xlsx", sheet = "input", skip = 1)

size_meta <- size_meta %>%
  rename(size = name)

modulus <- size_meta %>%
  left_join(raw_df, by = "size") %>%
  select(-sieve, -freq_hz) %>%
  filter(size != "Floccular") %>%
  pivot_longer(
    cols = c(G_1:G_3, G2_1:G2_3),
    names_to = c(".value", "replicate"),
    names_pattern = "(G2?)_(\\d)"
    ) 
  
modulus_subset <- modulus %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

# ------ Summarize mean across replicates for plotting ------

mod_summary_wide <- modulus %>%
  group_by(midpoint, size, freq_rad) %>%
  summarize(
    G_avg = mean(G),
    G_sd = sd(G),
    G2_avg = mean(G2),
    G2_sd = sd(G2),
    .groups = "drop"
  ) 

mod_summary_long <- mod_summary_wide %>%
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

mod_subset_sum_w <- mod_summary_wide %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

mod_subset_sum_l <- mod_summary_long %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

# ------ Correlation ------

# replicate level
res_storage <- cor.test(
  modulus_subset$G, 
  modulus_subset$midpoint, 
  method = "spearman", 
  exact = FALSE
  )

res_loss <- cor.test(
  modulus_subset$G2, 
  modulus_subset$midpoint, 
  method = "spearman",
  exact = FALSE
  )

# mean
res_storage_mean <- cor.test(
  mod_subset_sum_w$G_avg,
  mod_subset_sum_w$midpoint, 
  method = "spearman"
)

res_loss_mean <- cor.test(
  mod_subset_sum_w$G2_avg,
  mod_subset_sum_w$midpoint, 
  method = "spearman"
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
    legend.position = "right",
    strip.background = element_rect(
      colour = NA # facet label outline
      )
    )

p2 <- ggplot(mod_subset_sum_l, aes(x = size, y = avg, fill = measure)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  labs(
    title = "Frequency = 0.1 rad/s",
    x = "Size",
    y = "Modulus (kPa)"
  ) +
  scale_fill_manual(
    values = c("plum4", "lightgray")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

# arrange two plots into one column
p <- plot_grid(
  p1, p2,
  labels = "auto", ncol = 1, rel_widths = c(6.5, 5)
)

fname_out <- "./figures/Figure_1.tif"
ggsave(fname_out, plot = p, width = 6.5, height = 5, dpi = 300)

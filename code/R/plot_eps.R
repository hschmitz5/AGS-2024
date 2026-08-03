rm(list = ls())
library(tidyverse)
library(ggh4x)

# File names for concentration data
fname_pn    <- paste0("./data/EPS/PN_conc_ags.rds")
fname_polys <- paste0("./data/EPS/PS_conc_ags.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    filter(size != "Floccular") %>%
    group_by(extract, size) %>%
    summarize(
      avg = mean(C_VSS),
      sd = sd(C_VSS),
      .groups = "drop"
    ) %>%
    mutate(
      extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
      extract = factor(extract, levels = c("Tightly Bound", "Loosely Bound"))
    )
  }
# Apply function to each assay
PN <- group_data(fname_pn) 
PS <- group_data(fname_polys)

# Calculate PN + PS and PN/PS
df_wide <- left_join(
  PN %>% select(extract, size, PN_avg = avg, PN_sd = sd), 
  PS %>% select(extract, size, PS_avg = avg, PS_sd = sd), 
  by = c("extract", "size")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg,
    sd = NA
    ) 

# Combine into single data frame
df_conc <- bind_rows(
  'Protein (PN)' = PN,
  'Polysaccharide (PS)' = PS,
  'Total EPS (PN + PS)' = df_wide %>% select(extract, size, avg = total, sd),
  .id = "assay"
  ) %>%
  mutate(y_label = "\u00b5g/mgVSS") %>%  # Concentration
  select(y_label, assay, extract, size, avg, sd) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  mutate(
    y_label = "", # unitless
    assay = "PN/PS",
    sd = NA
  ) %>%
  select(y_label, assay, extract, size, avg = PNPS, sd) 

df_all <- bind_rows(df_conc, PNPS) %>%
  mutate(
    y_label = factor(y_label, levels = c("\u00b5g/mgVSS", "")),
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "Total EPS (PN + PS)", "PN/PS"))
  )


# ------ Plot ------

p <- ggplot(df_all, aes(x = size, y = avg, fill = assay)) +
  
  # Concentration Plots
  geom_col(
    data = subset(df_all, y_label == "\u00b5g/mgVSS"),
    position = "dodge",
    width = 0.8
  ) +
  geom_errorbar(
    data = subset(df_all, y_label == "\u00b5g/mgVSS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # PN/PS plots
  geom_col(
    data = subset(df_all, y_label == ""),
    width = 0.5
  ) +
  
  # Sizes
  ggh4x::facet_grid2(
    y_label ~ extract,
    scales = "free",
    switch = "y",
    independent = "y"
  ) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(),   
      scale_y_continuous(), 
      scale_y_continuous(breaks = c(0, 2, 4)), # PN/PS row
      scale_y_continuous(breaks = c(0, 2, 4))  # PN/PS row
    )
  ) +
  force_panelsizes(rows = c(1, 1/3), cols = c(1, 1)) +
  
  scale_fill_manual(
    values = c(
      "Polysaccharide (PS)" = "lightsalmon2",
      "Protein (PN)" = "lightblue",
      "Total EPS (PN + PS)" = "steelblue",
      "PN/PS" = "lightgray"
    )
  ) +
  
  labs(
    x = "Size",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank()
  )


fname_out <- "./figures/Figure_2.tif"
ggsave(fname_out, plot = p, width = 6.5, height = 3, dpi = 300)

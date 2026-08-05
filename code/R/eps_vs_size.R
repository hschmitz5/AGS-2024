rm(list = ls())
library(tidyverse)
library(ggh4x)

PN <- readRDS("./data/EPS/PN_conc_ags.rds") %>%
  select(extract, size, replicate, PN = C_VSS) 

PS <- readRDS("./data/EPS/PS_conc_ags.rds") %>%
  select(extract, size, replicate, PS = C_VSS)

eps <- left_join(PN, PS, by = c("extract", "size", "replicate")) %>%
  filter(size != "Floccular") %>%
  mutate(
    total = PN + PS, 
    ratio = PN/PS,
  )

summary_wide <- eps %>%
  group_by(extract, size) %>%
  summarize(
    # protein
    PN_avg = mean(PN), PN_sd = sd(PN),
    # polysaccharide
    PS_avg = mean(PS), PS_sd = sd(PS),
    # total 
    total_avg = mean(total), total_sd = sd(total),
    # PN/PS
    ratio_avg = mean(ratio), ratio_sd = sd(ratio),
    .groups = "drop"
  ) 

summary_long <- summary_wide %>%
  pivot_longer(
    cols = c(PN_avg, PN_sd,
             PS_avg, PS_sd,
             total_avg, total_sd,
             ratio_avg, ratio_sd),
    names_to = c("assay", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    y_label = if_else(assay == "ratio", "PN/PS", "\u00b5g/mgVSS"),
    # Write out LB and TB
    extract = factor(extract, levels = c("TB", "LB")),
    extract = recode(extract, "LB" = "Loosely Bound", "TB" = "Tightly Bound"),
    # write out PN, PS, etc
    assay = factor(assay, levels = c("PS", "PN", "total", "ratio")),
    assay = recode(assay, "PN" = "Protein (PN)", "PS" = "Polysaccharide (PS)",
                   "total" = "Total EPS (PN + PS)", "ratio" = "PN/PS")
    ) 


# ------ Plot ------

p <- ggplot(summary_long, aes(x = size, y = avg, fill = assay)) +
  
  # Concentration Plots
  geom_col(
    data = subset(summary_long, y_label == "\u00b5g/mgVSS"),
    position = "dodge",
    width = 0.8
  ) +
  geom_errorbar(
    data = subset(summary_long, y_label == "\u00b5g/mgVSS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # PN/PS plots
  geom_col(
    data = subset(summary_long, y_label == "PN/PS"),
    width = 0.5
  ) +
  geom_errorbar(
    data = subset(summary_long, y_label == "PN/PS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.5),
    width = 0.5/4
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

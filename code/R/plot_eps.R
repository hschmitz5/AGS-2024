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

# ------ Correlation ------

# define sample names
sz <- data.frame(
  size = c("S", "M", "L", "XL", "XXL"),
  midpoint = c(1.125, 1.7, 2.4, 3.4, 4.5)
)

df_wide <- df_wide %>%
  left_join(., sz, by = "size") %>%
  mutate(
    extract = recode(extract,"Loosely Bound" = "LB","Tightly Bound" = "TB")
  )

vars <- c("PN_avg", "PS_avg", "total", "PNPS") # rows

correlate_EPS <- function(df_wide, extract_type) {
  df_extract <- df_wide %>% filter(extract == extract_type)
  
  res_extract <- vars |>
    set_names() |>
    map_dfr(
      \(x) broom::tidy(
        cor.test(df_extract$midpoint, df_extract[[x]], method = "spearman")
      ),
      .id = "var"
    ) |>
    mutate(
      extract = extract_type,
      p.adj = p.adjust(p.value, method = "BH")
    ) |>
    select(extract, var, p.value, p.adj, estimate)
}

res_LB <- correlate_EPS(df_wide, "LB")
res_TB <- correlate_EPS(df_wide, "TB")

res <- rbind(res_TB, res_LB)

# Differential Abundance comparison
# Convert EPS concentrations to differential abundance (relative to S)
DA_eps <- df_wide %>%
  filter(extract == "LB") %>%
  mutate(
    PN_diff = PN_avg - PN_avg[size == "S"][1],
    PS_diff = PS_avg - PS_avg[size == "S"][1],
    total_diff = total - total[size == "S"][1],
    PNPS_diff = PNPS - PNPS[size == "S"][1]
  ) %>%
  filter(size != "S") %>%
  select(extract, size, PN_diff, PS_diff, total_diff, PNPS_diff)

# Community Correlation
DA_comm <- readRDS("./data/DA/DA_genus_processed.rds") %>%
  filter(Genus == "Ca_Contendobacter")

vars <- c("PN_diff", "PS_diff", "total_diff", "PNPS_diff")
DA_test <- vars |>
  set_names() |>
  map_dfr(
    \(x) broom::tidy(
      cor.test(DA_comm$lfc, DA_eps[[x]], method = "spearman")
    ),
    .id = "var"
  ) |>
  mutate(
    p.adj = p.adjust(p.value, method = "BH")
  ) |>
  select(var, p.value, p.adj, estimate)

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


fname_out <- "./figures/EPS.png"
ggsave(fname_out, plot = p, width = 6.5, height = 3, dpi = 300)

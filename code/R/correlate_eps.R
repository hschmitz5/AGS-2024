rm(list = ls())
library(readxl)
library(tidyverse)
library(patchwork)
source("./code/R/01_load_ps.R")
rm(ps) # do not confuse ps with PS

# ------ Rheometry data ------

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
    size = factor(size, levels = size_meta$name),
    measure = factor(measure, levels = c("G", "G2"))
  )

modulus_subset <- modulus %>%
  filter(
    measure == "G", # G2
    freq_rad == 0.1
    ) %>%
  select(size, avg) 

# ------ EPS data ------

# File names for concentration data
fname_pn    <- paste0("./data/EPS/PN_conc_ags.rds")
fname_polys <- paste0("./data/EPS/PS_conc_ags.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    filter(size != "Floccular") %>%
    group_by(size, extract) %>%
    summarize(
      avg = mean(C_VSS),
      sd = sd(C_VSS),
      .groups = "drop"
    )
}
# Apply function to each assay
PN <- group_data(fname_pn) 
PS <- group_data(fname_polys)

# Change name for join
size_meta <- size_meta %>%
  rename(size = name) %>%
  select(-ranges)

# Calculate PN/PS
df_wide <- left_join(
  PN %>% select(size, extract, PN_avg = avg), 
  PS %>% select(size, extract, PS_avg = avg), 
  by = c("size", "extract")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg
  ) %>%
  left_join(., size_meta, by = "size") %>%       # join midpoint
  left_join(., modulus_subset, by = "size") %>%  # join modulus
  rename(modulus = avg) 

# ------ Differential Abundance data ------

# Community Correlation
DA_comm <- readRDS("./data/DA/DA_genus_processed.rds") %>%
  filter(Genus == "Ca_Contendobacter")

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

# ------ Correlate Size ------

var2 <- c("PN_avg", "PS_avg", "total", "PNPS") # rows

correlate_eps <- function(df_wide, var1, extract_type) {

  df_extract <- df_wide %>% filter(extract == extract_type)

  res_extract <- var2 |>
    set_names() |>
    map_dfr(
      \(x) broom::tidy(
        cor.test(df_extract[[var1]], df_extract[[x]], method = "spearman")
      ),
      .id = "var"
    ) |>
    mutate(
      extract = extract_type,
      p.adj = p.adjust(p.value, method = "BH")
    ) |>
    select(extract, var, p.value, p.adj, estimate)

}

res_midpoint <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(df_wide, "midpoint", .x))

res_modulus <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(df_wide, "modulus", .x))



### Differential Abundance

vars <- c("PN_diff", "PS_diff", "total_diff", "PNPS_diff")

res_DA <- vars |>
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

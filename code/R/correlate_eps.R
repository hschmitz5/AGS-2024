rm(list = ls())
library(readxl)
library(tidyverse)
source("./code/R/01_load_ps.R")
rm(ps) # do not confuse ps with PS

# ------ Rheometry data ------

fname_in <- "./data/Rheometry_Nov_2024.xlsx"

mod_summary <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  select(-sieve, -freq_hz) %>%
  filter(size != "Floccular") %>%
  pivot_longer(
    cols = c(G_1:G_3, G2_1:G2_3),
    names_to = c(".value", "replicate"),
    names_pattern = "(G2?)_(\\d)"
  ) %>%
  # Summarize means
  group_by(size, freq_rad) %>%
  summarize(
    G_avg = mean(G),
    G_sd = sd(G),
    G2_avg = mean(G2),
    G2_sd = sd(G2),
    .groups = "drop"
  ) %>%
  mutate(
    size = factor(size, levels = size_meta$name)
  )

mod_summary_subset <- mod_summary %>%
  filter(
    freq_rad == 0.1
    ) %>%
  select(size, G_avg, G2_avg) 

# ------ EPS data ------

# Change name for join
size_meta <- size_meta %>%
  rename(size = name) 

PN <- readRDS("./data/EPS/PN_conc_ags.rds") %>%
  select(extract, size, replicate, PN = C_VSS) %>%
  filter()

PS <- readRDS("./data/EPS/PS_conc_ags.rds") %>%
  select(extract, size, replicate, PS = C_VSS)

eps_conc <- left_join(PN, PS, by = c("extract", "size", "replicate")) %>%
  filter(size != "Floccular") %>%
  # join midpoint
  left_join(., size_meta, by = "size") %>%
  select(extract, midpoint, size, replicate, PN, PS) %>%
  # Calculate total and ratio
  mutate(
    total = PN + PS, 
    ratio = PN/PS,
    size = factor(size, levels = size_meta$size)
  ) 

eps_summary <- eps_conc %>%
  group_by(extract, midpoint, size) %>%
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
  ) %>%
  # join modulus
  left_join(., mod_summary_subset, by = "size") 
  

# ------ ANOVA ------

df <- eps_conc %>%
  filter(extract == "LB")

# ANOVA
mod_PN    <- aov(PN ~ size, data = df)

summary(mod_PN)

bartlett.test(PN ~ size, data = df)

# par(mfrow = c(2, 2))
# plot(mod)
  
# ------ Correlate Size ------

correlate_eps <- function(df_wide, var1, extract_type) {

  df_extract <- df_wide %>% filter(extract == extract_type)

  res_extract <- var2 |>
    set_names() |>
    map_dfr(
      \(x) broom::tidy(
        cor.test(
          df_extract[[var1]], 
          df_extract[[x]], 
          method = "spearman",
          exact = FALSE
          )
      ),
      .id = "var"
    ) |>
    mutate(
      extract = extract_type,
      p.adj = p.adjust(p.value, method = "BH")
    ) |>
    select(extract, var, p.value, p.adj, estimate)

}

# replicate level variables
var2 <- c("PN", "PS", "total", "ratio")

# EPS vs size (replicate level)
res_midpoint_rep <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(eps_conc, "midpoint", .x))

# mean level variables
var2 <- c("PN_avg", "PS_avg", "total_avg", "ratio_avg")

# EPS vs size (mean)
res_midpoint <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(eps_summary, "midpoint", .x))

# Mean EPS vs mean modulus
res_modulus_G <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(eps_summary, "G_avg", .x))

res_modulus_G2 <- c("TB", "LB") |>
  map_dfr(~ correlate_eps(eps_summary, "G2_avg", .x))


# ------ Differential Abundance data ------

# Community Correlation
DA_comm <- readRDS("./data/DA/DA_genus_processed.rds") %>%
  filter(Genus == "Ca_Contendobacter")

# Convert EPS concentrations to differential abundance (relative to S)
DA_eps <- eps_summary %>%
  filter(extract == "LB") %>%
  mutate(
    PN_diff    = PN_avg - PN_avg[size == "S"][1],
    PS_diff    = PS_avg - PS_avg[size == "S"][1],
    total_diff = total_avg - total_avg[size == "S"][1],
    ratio_diff = ratio_avg - ratio_avg[size == "S"][1]
  ) %>%
  filter(size != "S") %>%
  select(extract, size, PN_diff, PS_diff, total_diff, ratio_diff)


### Differential Abundance

vars <- c("PN_diff", "PS_diff", "total_diff", "ratio_diff")

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

rm(list = ls())
library(readxl)
library(tidyverse)
library(patchwork)
source("./code/R/01_load_ps.R")
rm(ps) # do not confuse ps with PS

# ------ Rheometry data ------

# define sample names
sz <- data.frame(
  name = c("S", "M", "L", "XL", "XXL")
)

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
    size = factor(size, levels = sz$name),
    measure = factor(measure, levels = c("G", "G2"))
  )

modulus_subset <- modulus %>%
  filter(
    measure == "G2",
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

# Calculate PN/PS
df_wide <- left_join(
  PN %>% select(size, extract, PN_avg = avg), 
  PS %>% select(size, extract, PS_avg = avg), 
  by = c("size", "extract")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg
  )

# ------ Correlation ------

vars <- c("PN_avg", "PS_avg", "total", "PNPS") # rows

correlate_EPS <- function(df_wide, extract_type) {
  df_extract <- df_wide %>% filter(extract == extract_type)
  
  res_extract <- vars |>
    set_names() |>
    map_dfr(
      \(x) broom::tidy(
        cor.test(modulus_subset$avg, df_extract[[x]], method = "spearman")
        ),
      .id = "var"
    ) |>
    mutate(
      extract = extract_type,
      p.adj = p.adjust(p.value, method = "BH")
    ) |>
    select(extract, var, p.value, p.adj, estimate)
}

res_TB <- correlate_EPS(df_wide, "TB")
res_LB <- correlate_EPS(df_wide, "LB")

res <- rbind(res_TB, res_LB)

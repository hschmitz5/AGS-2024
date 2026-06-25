# NOTE: model forces intercept to intersect (0,0)

# Verify quality of predicted output 
# by verifying that the data is within the working range,
# and check the R-squared in the model summary

rm(list = ls())
library(readxl)
library(writexl)
library(tidyverse)
library(ggplot2)

# change to process each sheet
sheet_name <- "PS"  # protein (PN) or polysaccharide (PS)

# File name for absorbance data
fname_in   <- "./data/EPS/EPS_absorbance.xlsx"

# Used to convert concentration to include TSS/VSS
extract_volume = 10 # mL

# define sample names
sz <- data.frame(
  sieve = c(40, 20, 14, 10, 7, 5),
  name  = c("XS", "S", "M", "L", "XL", "XXL")
)

# ----------------------------------------


# Read in absorbance data and rename size
df <- read_excel(fname_in, sheet = sheet_name, skip = 1) %>%
  left_join(sz, by = c("size" = "sieve")) %>%
  mutate(size = factor(name, levels = sz$name)) %>%
  select(-name)

# define sample data (excludes standards)
sam <- df %>%
  filter(is.na(C)) %>%
  select(-C) %>%
  # initialize data 
  mutate(
    C0 = NA_real_,
    C_VSS = NA_real_) 

# set poly_degree based on sheet_name
poly_degree <- switch(sheet_name,
                      "PN" = 3,
                      "PS" = 2,
                      stop("invalid sheet_name"))

# loop over each dataset
for (d_set in unique(df$dataset)) {
  # define sample index corresponding to dataset
  idx <- sam$dataset == d_set
  # define standards for dataset
  std <- df %>%
    filter(
      dataset == d_set,
      is.finite(C)
    ) %>%
    select(sample,A,C)
  
  ## fit concentration as a function of absorbance
  # 0 forces the equation to intercept (0,0)
  model <- lm(C ~ 0 + poly(A, poly_degree, raw = TRUE), data = std)
  msum <- summary(model)
  
  # predict concentrations
  sam$C0[idx] <- predict(model, newdata = sam[idx,])
  # convert units to VSS
  sam$C_VSS[idx] <- sam$C0[idx]*extract_volume/sam$VSS[idx]
  
  #### Plot Fit Data
  
  # fit line for plotting
  A_seq <- seq(min(std$A), max(std$A), length.out = 200)
  fit_df <- data.frame(
    A = A_seq,
    C = predict(model, newdata = data.frame(A = A_seq))
  )
  
  # plot
  ggplot() +
    geom_line(data = fit_df, aes(C, A, color = "Fit"), linewidth = 1) +
    geom_point(data = std, aes(C, A, color = "Standards")) +
    geom_point(data = sam %>% filter(dataset == d_set), 
               aes(C0, A, color = "Samples"), shape = 3) +
    scale_color_manual(
      values = c(
        "Fit" = "lightblue",
        "Standards" = "black",
        "Samples" = "red"
      ),
      labels = c(
        "Fit" = paste0("Fit (R² = ", round(msum$r.squared, 4), ")"),
        "Standards" = "Standards",
        "Samples" = "Samples"
      )
    ) +
    labs(
      color = NULL,
      x = bquote("Concentration [" * mu * "g" * .(sheet_name) * "/mL]"),
      y = "Absorbance"
    ) +
    theme_minimal(base_size = 12) +
    theme(aspect.ratio = 0.7)
  
  fit_plot <- paste0("./figures/EPS/",sheet_name,"_fit_set_",d_set,".png")
  ggsave(fit_plot, height = 2.5, width = 6, dpi = 600)
}

sam_output <- sam %>%
  dplyr::select(size, replicate, extract, C_VSS)

# save the sample data
saveRDS(sam_output, file = paste0("./data/EPS/",sheet_name,"_conc.rds"))

sam_sort <- sam_output %>%
  arrange(extract)

write_xlsx(sam_sort, path = paste0("./data/EPS/",sheet_name,"_conc.xlsx"))

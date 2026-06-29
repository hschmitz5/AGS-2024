rm(list = ls())
library(tidyverse)
library(patchwork)

# File names for concentration data
fname_pn    <- paste0("./data/EPS/PN_conc.rds")
fname_polys <- paste0("./data/EPS/PS_conc.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    filter(size != "XS") %>%
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

# Calculate PN + PS and PN/PS
df_wide <- left_join(
  PN %>% select(size, extract, PN_avg = avg, PN_sd = sd), 
  PS %>% select(size, extract, PS_avg = avg, PS_sd = sd), 
  by = c("size", "extract")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg,
    sd = NA
    ) 

# Combine into single data frame
df <- bind_rows(
  'Protein (PN)' = PN,
  'Polysaccharide (PS)' = PS,
  'PN + PS' = df_wide %>% select(size, extract, avg = total, sd),
  .id = "assay"
  ) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "PN + PS"))
    ) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  select(size, extract, avg = PNPS) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound")
    ) 

# Determine maximum avg + sd
max_y1 <- df %>%
  filter(assay != "PN + PS") %>%
  summarise(max_y = max(avg + sd)) %>%
  pull(max_y)

max_y2 <- df %>%
  filter(assay == "PN + PS") %>%
  summarise(max_y = max(avg)) %>%
  pull(max_y)

max_y <- ceiling(
  max(max_y1, max_y2)
  )

# -------------------------------

# Make Plot

p <- ggplot(data = df, aes(x = size, y = avg, fill = assay)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~extract) + 
  ylim(0, max_y) +
  labs(
    y = expression(paste(mu, "g/mgVSS")),
    x = NULL,
    fill = NULL # legend titles
  ) +
  scale_fill_manual(
    values = c(
      "Polysaccharide (PS)" = "lightsalmon2",
      "Protein (PN)"        = "lightblue",
      "PN + PS"      = "steelblue" 
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      colour = NA # facet label outline
    ),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) 

annot <- ggplot(data = PNPS, aes(x = size, y = avg, fill = "PN/PS")) +
  geom_col(position = "dodge", width = 0.5) +
  facet_wrap(~extract) +
  scale_y_continuous(
    breaks = c(0, 2, 4)
    ) +
  labs(
    x = "Size",
    y = NULL, 
    fill = NULL
    ) +
  scale_fill_manual(
    values = "lightgray",
    labels = expression(frac(PN, PS))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text  = element_blank(),
    legend.justification = "left"
  )

p2 <- p / annot +
  plot_layout(heights = c(4, 1.5))

fname_out <- "./figures/EPS.png"
ggsave(fname_out, plot = p2, width = 6.5, height = 3, dpi = 300)

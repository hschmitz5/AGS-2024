rm(list = ls())
library(tidyverse)
library(patchwork)

# File names for concentration data
fname_pn    <- paste0("./data/EPS/PN_conc.rds")
fname_polys <- paste0("./data/EPS/PS_conc.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
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

# Combine into single data frame
df <- bind_rows(
  Protein = PN,
  Polysaccharide = PS,
  .id = "assay"
  ) %>%
  group_by(size, extract) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
    assay = factor(assay, levels = c("Protein", "Polysaccharide"))
    ) 

# Calculate PN/PS
df_wide <- PN %>% 
  rename(PN_avg = avg, PN_sd = sd) %>%
  left_join(
    PS %>% rename(PS_avg = avg, PS_sd  = sd),
    by = c("size", "extract")
  ) %>%
  mutate(
    PNPS = PN_avg/PS_avg,
    total_avg = PN_avg + PS_avg,
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound")
  )

tot <- df_wide %>%
  transmute(assay = "Total EPS", size, extract, avg = total_avg) 

# Determine maximum avg + sd
max_y <- ceiling(
  max(df$avg + df$sd)
)

# -------------------------------

# Make Plot

p <- ggplot(data = df, aes(x = size, y = avg, fill = assay)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.6),
    width = 0.2
  ) +
  geom_point(data = tot, aes(x = size, y = avg, color = assay), shape = 20) +
  facet_wrap(~extract, nrow=1, strip.position = "top") +
  ylim(0, max_y) +
  labs(
    x = "Size",
    y = expression(paste(mu, "g/mgTSS")),
    color = NULL,
    fill = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Protein" = "gray",
      "Polysaccharide" = "lightblue"
    )
  ) +
  scale_color_manual(
    values = c("Total EPS" = "black")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(
      fill = "white",      # facet label fill
      colour = "lightgray" # facet label outline
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) 
  
annot <- ggplot(data = df_wide) +
  geom_tile(aes(x = size, y = "PN/PS", fill = round(PNPS,1))) +
  geom_text(aes(x = size, y = "PN/PS", label = round(PNPS,1))) +
  scale_fill_gradient(low="white", high="lightgray") +
  facet_wrap(~extract, nrow=1, strip.position = "top") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid  = element_blank(),
    axis.title  = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks  = element_blank(),
    strip.text  = element_blank()
  )

p2 <- p / annot +
  plot_layout(heights = c(4, 1)) 

fname_out <- "./figures/EPS.png"
ggsave(fname_out, plot = p2, width = 8, height = 4, dpi = 600)

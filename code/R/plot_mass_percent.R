source("./code/R/01_load_ps.R")

# Data -------------------------------------------------------------

labels <- c("< 0.21", "0.21 - 0.43", "0.43 - 0.60", "0.60 - 1.4", 
            "1.4 - 2.0", "2.0 - 2.8", "2.8 - 4.0", "> 4.0")

x_edges <- c(0, 0.21, 0.43, 0.60, 1.4, 2.0, 2.8, 4.0, 5.0)   # sieve boundaries

mass_percent <- c(27.29, 19.05, 5.05, 11.88, 7.59, 9.85, 10.62, 8.67)

colors <- c("#606060", "white", "gray", met.brewer(size_pal, n_sizes))

df <- data.frame(
  xmin = x_edges[-length(x_edges)],  # removes last element
  xmax = x_edges[-1],                # removes first element 
  ymin = 0,
  ymax = mass_percent,
  label = factor(labels, levels = labels),
  color = colors 
)

# Plot --------------------------------------------------------------

p <- ggplot(df) +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label), 
    color = "black", linewidth = 0.2  # outline
  ) +
  scale_fill_manual(values = colors) +
  # control x-axis tick labels
  scale_x_continuous(
    breaks = x_edges[-length(x_edges)],
    labels = x_edges[-length(x_edges)]
  ) +
  labs(
    x = "Biomass Diameter [mm]",
    y = "Mass Percentage [%]",
    fill = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor.x = element_blank()
  ) +
  guides(
    fill = guide_legend(byrow = TRUE)
  )

# Save figure -------------------------------------------------------

fname <- "./figures/mass-percent.png"
ggsave(fname, p, width = 6, height = 5, dpi = 300)

# Display
print(p)

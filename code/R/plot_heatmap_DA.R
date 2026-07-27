rm(list = ls())
source("./code/R/01_load_ps.R")
library(ComplexHeatmap)
library(circlize) # for colorRamp2

n_display_low <- 30

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.6 # same as cell_h

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

rel_ab_cutoff <- 0.5
p_threshold <- 0.05

# number of rows to show
n_show <- 30

# ------ Process Data -----

data_df <- readRDS("./data/DA/DA_genus_processed.rds") 

# Genera mapped in Relative Abundance
rel_names <- get_rel_wide(ps) %>%
  # Arrange taxa from largest to smallest abundance
  mutate(row_sum = rowSums(.)) %>%
  arrange(desc(row_sum)) %>%
  # Keep the top n_show
  head(., n = n_show) %>%
  rownames()

# Names of genera with metabolism annotation
# m_names <- get_metabolism(data_df) %>%
#   filter(if_any(everything(), ~ !is.na(.))) %>%
#   rownames()

data_mat <- data_df %>%
  filter(Genus %in% rel_names) %>%
  dplyr::select(-robust) %>%
  pivot_wider(names_from = size, values_from = lfc) %>%
  column_to_rownames(var = "Genus") %>%
  as.matrix()

# Metabolism
m_df <- as.data.frame(data_mat) %>%
  rownames_to_column(var = "Genus") %>%
  get_metabolism() 

# ---- Plotting

# metabolism annotation
m_colors  <- c("P" = "#66C24A", "V" = "#EAEC3F") 
m_annot <- rowAnnotation(
  df = m_df,
  gap = unit(2, "mm"),
  # column names
  annotation_name_side = "top",
  annotation_name_rot = 60,
  # color
  col = col_list <- setNames(
    rep(list(m_colors), ncol(m_df)),
    colnames(m_df)
  ),
  na_col = NA, # no color for NA
  # legend
  show_legend = FALSE
)
# metabolism legend
metab_lgd <- Legend(
  title = "Functional Group",
  labels = c("Positive", "Variable"),
  legend_gp = gpar(fill = m_colors),
  nrow = 1,
  title_position = "leftcenter",
  title_gap = unit(3, "mm")
)

# Dimensions
n_cols <- ncol(data_mat)
n_rows <- nrow(data_mat)

# Labels
row_labels <- rownames(data_mat)
italic_rows <- !grepl("^(Unk|midas)", row_labels) 
row_fontface <- ifelse(italic_rows, "italic", "plain")

# Heatmap
ht_colors <- colorRamp2(
  c(min(data_mat), 0, max(data_mat)), 
  c("dodgerblue4", "white", "red3")
)

breaks_display <- c(-3, -1.5, 0, 1.5, 3)

heatmap_lgd <- Legend(
  col_fun = ht_colors,
  labels = breaks_display,
  at = breaks_display,
  title = "Log Fold-Change (Relative to S)",
  direction = "horizontal",
  legend_width = unit(5.9, "cm")
)

ht <- Heatmap(
  data_mat,
  # columns
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  cluster_columns = FALSE, # changes sample order
  # heatmap legend
  col = ht_colors, 
  show_heatmap_legend = FALSE, 
  # Annotations
  right_annotation = m_annot,
  # Display size
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize, fontface = "bold")
)

pd = packLegend(heatmap_lgd, metab_lgd, direction = "horizontal", 
                max_width = unit(10, "cm"), row_gap = unit(5, "mm"))

fname  <- "./figures/genus_level_DA.png"
# Draw combined heatmap
png(fname,
    width = 5.5,  # width in inches; can adjust
    height = 8.3, # height in inches; can adjust
    units = "in", res = 300)
draw(ht) 
draw(pd, x = unit(0.1, "npc"), y = unit(0.03, "npc"), just = c("left", "bottom"))
dev.off()

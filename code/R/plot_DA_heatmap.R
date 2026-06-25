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
  # column names
  annotation_name_side = "bottom",
  annotation_name_rot = -60,
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
lgd <- Legend(
  title = NULL,
  labels = c("Positive", "Variable"),
  legend_gp = gpar(fill = m_colors),
  nrow = 2,
  row_gap = unit(3, "mm")
)

# Dimensions
n_cols <- ncol(data_mat)
n_rows <- nrow(data_mat)

# Labels
row_labels <- rownames(data_mat)
italic_rows <- !grepl("^(Unk|midas)", row_labels) 
row_fontface <- ifelse(italic_rows, "italic", "plain")

# Heatmap
col_fun <- colorRamp2(
  c(min(data_mat), 0, max(data_mat)), 
  c("dodgerblue4", "white", "red3")
)

breaks_display <- c(-3, -1.5, 0, 1.5, 3)

ht <- Heatmap(
  data_mat,
  # columns
  column_title = NULL, #"Differential Abundance (log fold change)", 
  cluster_columns = FALSE, # changes sample order
  column_names_rot = 0,
  column_names_centered = TRUE,
  # heatmap legend
  col = col_fun, 
  heatmap_legend_param = list(
    at = breaks_display,
    labels = breaks_display,
    title = NULL, #"log fold change", 
    direction = "horizontal",
    legend_width = unit(4, "cm")
  ),
  # # Annotations
  right_annotation = m_annot,
  # Display size
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize)
)

fname  <- "./figures/DA_heatmap.png"
# Draw combined heatmap
png(fname,
    width = 7,  # width in inches; can adjust
    height = 6.5, # height in inches; can adjust
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "top") #, annotation_legend_side = "top") 
draw(lgd, x = unit(0.66, "npc"), y = unit(0.98, "npc"), just = c("right", "top"))
dev.off()

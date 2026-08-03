rm(list = ls())
library(ComplexHeatmap)
library(circlize) # for colorRamp2
source("./code/R/01_load_ps.R")

# number of rows to show
n_show <- 30 

write2excel <- 0

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.2 

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

# ------ Define Data ------

# Relative Abundance
rel_df <- get_rel_wide(ps) %>%
  # Arrange taxa from largest to smallest abundance
  mutate(row_sum = rowSums(.)) %>%
  arrange(desc(row_sum)) %>%
  # Keep the top n_show
  head(., n = n_show) %>%
  dplyr::select(-row_sum) 

# sum of relative abundance in each sample
rel_sum <- colSums(rel_df)

# Convert to log
log_mat <- as.matrix(rel_df) %>%
  { log10(.) }

# Metabolism
m_df <- rel_df %>%
  rownames_to_column(var = "Genus") %>%
  get_metabolism() 

# DA Taxa
DA_genera <- readRDS("./data/DA/DA_genus_processed.rds") %>%
  pull(Genus) %>%
  unique()

#### ---- Plotting ------

# size annotation
size_annot <- HeatmapAnnotation(
  sz = anno_block(
    gp = gpar(
      fill = c(rep("lightgray", n_sizes)), # all the same color
      col = NA # removes border
    ),
    labels = size_meta$name,
    labels_gp = gpar(
      col = c(rep("black", n_sizes)), 
      fontsize = col_fontsize
    )
  )
)

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

## Relative Abundance

# Dimensions
n_cols <- ncol(log_mat)
n_rows <- nrow(log_mat)
split = rep(1:n_sizes, each = n_replicates)

# Labels
row_labels <- rownames(log_mat)
# italicize classified genera
italic_rows <- !grepl("^(Unk|midas)", row_labels) 
# bold DA taxa
bold_rows <- row_labels %in% DA_genera
# Apply
row_fontface <- ifelse(
  bold_rows & italic_rows, "bold.italic",
  ifelse(bold_rows, "bold",
         ifelse(italic_rows, "italic", "plain"))
)

# Relative Abundance legend 
break_values <- c(round(min(rel_df),2), 1, 10, round(ceiling(max(rel_df))))
breaks_log_display <- log10(break_values) # Log (%)
# Add % symbol to largest break
breaks_rel_display <- replace(
  as.character(break_values),
  length(break_values),
  paste0(tail(break_values, 1), "%")
)
# Legend Colors
ht_colors <- colorRamp2(
  breaks <- seq(
    log10(break_values[1]),
    log10(tail(break_values, 1)),
    length.out = 10
  ), 
  met.brewer(taxa_pal, 10)
)

heatmap_lgd <- Legend(
  col_fun = ht_colors,
  labels = breaks_rel_display,
  at = breaks_log_display,
  title = "Relative Abundance",
  direction = "horizontal",
  legend_width = unit(7.2, "cm")
)

ht <- Heatmap(
  log_mat,
  # columns
  cluster_columns = FALSE, # changes sample order
  show_column_names = FALSE,
  column_split = split, # put a gap between sizes
  column_title = NULL,  # sample names
  # heatmap legend
  col = ht_colors, 
  show_heatmap_legend = FALSE, 
  # Annotations
  top_annotation = size_annot,
  right_annotation = m_annot, 
  # Display size
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize)
)

pd = packLegend(heatmap_lgd, metab_lgd, direction = "horizontal", 
                max_width = unit(10, "cm"), row_gap = unit(5, "mm"))

# Figure output location
fname_rel <- "./figures/heatmap_rel_ab.tif"

# Draw combined heatmap
png(fname_rel,
    width = 6.2,  # width in inches; can adjust
    height = 9.5, # height in inches; can adjust
    units = "in", res = 300)
draw(ht)
draw(pd, x = unit(0.1, "npc"), y = unit(0.03, "npc"), just = c("left", "bottom"))
dev.off()

## Check what percent of relative abundance is included in plot
message(paste("Heatmap Min:", round(min(rel_sum), 2), "%"))
message(paste("Heatmap Max:", round(max(rel_sum), 2), "%"))

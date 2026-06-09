rm(list = ls())
library(ComplexHeatmap)
source("./code/R/01_load_ps.R")
source("./code/R/02_metab_and_DA.R")

# Figure output location
fname_rel <- "./figures/rel_ab_heatmap.png"

# number of rows to show
n_show <- 30

write2excel <- 0

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.2 

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

# Pseudo count to add when converting to log
# (choose based on detection limit)
pseudo <- 1e-6  

# ------ Define Data ------

# Relative Abundance
rel_wide <- get_rel_wide(ps) %>%
  # Arrange taxa from largest to smallest abundance
  mutate(row_sum = rowSums(.)) %>%
  arrange(desc(row_sum)) %>%
  # Keep the top n_show
  head(., n = n_show) %>%
  dplyr::select(-row_sum) 

## Check what percent of relative abundance is included in plot
# sum per sample
rel_sum <- colSums(rel_wide)

# Convert to log
data_mat <- as.matrix(rel_wide) %>%
  { log10(. + pseudo) } 

# Metabolism
m_df <- rel_wide %>%
  rownames_to_column(var = "Genus") %>%
  get_metabolism(., metab_fname) 

# DA taxa
DA_taxa <- get_ancom_taxa(ancom_fname, ps, p_threshold)


#### ---- Plotting ------

n_cols <- ncol(data_mat)
n_rows <- nrow(data_mat)
split = rep(1:n_sizes, each = n_replicates)

# Labels
row_labels <- rownames(data_mat)
# italicize species + _g, but not _f/_o/_c/_p
italic_rows <- !grepl("^(Unk|midas)", row_labels) 
# bold significant taxa
bold_rows <- row_labels %in% DA_taxa
# Apply
row_fontface <- ifelse(
  bold_rows & italic_rows, "bold.italic",
  ifelse(bold_rows, "bold",
         ifelse(italic_rows, "italic", "plain"))
)

## Colors
# Relative Abundance
ht_colors <- met.brewer(taxa_pal, type = "continuous")
# Size
sz_colors <- c(rep("lightgray", n_sizes)) 
# Metabolism
m_colors  <- c("P" = "#66C24A", "V" = "#EAEC3F") 

# size annotation
size_annot <- HeatmapAnnotation(
  sz = anno_block(
    gp = gpar(
      fill = sz_colors,
      col = NA # removes border
    ),
    labels = size$name,
    labels_gp = gpar(
      col = c(rep("black", n_sizes)), 
      fontsize = col_fontsize
    )
  )
)

# metabolism annotation
m_annot <- rowAnnotation(
  df = m_df,
  # color
  col = col_list <- setNames(
    rep(list(m_colors), ncol(m_df)),
    colnames(m_df)
  ),
  na_col = NA, # no color for NA
  # legend
  show_legend = c(TRUE, c(rep(FALSE, length(m_df)-1))),  # Only first column contributes
  annotation_legend_param = list(
    title = NULL, 
    at = c("P", "V"),
    labels = c("Positive", "Variable"),
    nrow = 2
  ),
  annotation_name_side = "bottom",
  annotation_name_rot = -60
)

ht <- Heatmap(
  data_mat,
  # columns
  column_title = NULL, # "Relative Abundance",
  cluster_columns = FALSE, # changes sample order
  show_column_names = FALSE,
  column_split = split, # put a gap between sizes
  # heatmap legend
  show_heatmap_legend = TRUE, 
  col = ht_colors,
  heatmap_legend_param = list(
    title = "Log (%)", 
    direction = "horizontal",
    title_position = "lefttop"
  ),
  # Annotations
  bottom_annotation = size_annot,
  right_annotation = m_annot, 
  # Display size
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize)
)

# Draw combined heatmap
png(fname_rel,
    width = 7,  # width in inches; can adjust
    height = 8, # height in inches; can adjust
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top") 
dev.off()

# min and max of all sample sums
message(paste("Heatmap Min:", round(min(rel_sum), 2), "%"))
message(paste("Heatmap Max:", round(max(rel_sum), 2), "%"))

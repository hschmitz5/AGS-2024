rm(list = ls())
library(ComplexHeatmap)
source("./code/R/01_load_ps.R")
source("./code/R/02_metab_and_DA.R")

# Figure output location
fname_rel <- "./figures/rel_ab_heatmap.png"

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

# --------- Data ------------

#### Rename: agglomerate names when multiple ASVs are differentially abundant or not
rel_wide <- get_rel_agglom(ps, ancom_fname, rel_ab_cutoff, p_threshold)

# Convert to log
data_mat <- rel_wide %>%
  select(-Genus, -DA) %>%
  column_to_rownames("OTU") %>%
  as.matrix() %>%
  { log10(. + pseudo) }

# Names of DA taxa
DA_taxa_renamed <- rel_wide %>%
  filter(DA == "T") %>%
  pull(OTU)

#### Load metabolism
m_df <- as.data.frame(get_metabolism(rel_wide, metab_fname)) %>%
  rename(`NO2 reduction` = `Nitrite reduction`)

# ---- Plotting
n_cols <- ncol(data_mat)
n_rows <- nrow(data_mat)
split = rep(1:n_sizes, each = n_replicates)

# Labels
row_labels <- rownames(data_mat)
# italicize species + _g, but not _f/_o/_c/_p
italic_rows <- grepl("_(g|s)(?:_|$|-)", row_labels)
# bold significant taxa
bold_rows <- row_labels %in% DA_taxa_renamed
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
  column_title = "Relative Abundance",
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
    width = 8,  # width in inches; can adjust
    height = 10.5, # height in inches; can adjust
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top") 
dev.off()


#### Check what percent of relative abundance is included in plot
# sum per sample
rel_sum <- rel_wide %>%
  dplyr::select(-Genus, -DA) %>%
  column_to_rownames(var = "OTU") %>%
  colSums() 
# min and max of all sample sums
message(paste("Heatmap Min:", round(min(rel_sum), 2), "%"))
message(paste("Heatmap Max:", round(max(rel_sum), 2), "%"))


#### Export to Excel

if (write2excel == 1) {
  new_m <- m_df %>%
    # tf is true if any values in row are defined
    mutate(tf = as.integer(if_any(everything(), ~ !is.na(.x)))) %>%
    rownames_to_column(var = "OTU")
  
  full_df <- left_join(rel_wide, new_m, by = "OTU") %>%
    filter(tf == 1) %>%
    dplyr::select(-DA, -tf) %>%
    relocate(where(is.numeric), .after = where(is.character)) %>%
    arrange(Genus)
  
  library(writexl)
  write_xlsx(full_df, path = "./data/rel_ab_metab_high_ab.xlsx")
}
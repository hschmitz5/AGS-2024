rm(list = ls())
source("./code/R/01_load_ps.R")
source("./code/R/02_metab_and_DA.R")
library(writexl)
library(ComplexHeatmap)
library(circlize) # for colorRamp2

write2excel <- FALSE
n_display_low <- 30

fname_excel <- "./data/ANCOM_ASV.xlsx"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.6 # same as cell_h

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

rel_ab_cutoff <- 0.5

# ------ Process Data -----

taxonomy <- get_taxonomy(ps)

# define taxa in which at least one sample has abundance > rel_ab_cutoff
high_ab_taxa <- get_rel(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

# load differential abundance data
output <- readRDS(ancom_fname)

# At least one samples must be significant and pass sensitivity analysis
DA_taxa <- output$res %>%
  rename(OTU = taxon) %>%
  filter(OTU %in% high_ab_taxa) %>%
  # Combines diff_size* and passed_ss* together
  pivot_longer(
    cols = matches("q_size\\.name|passed_ss_size\\.name"),
    names_to = c(".value","size"),
    names_pattern = "(q|passed_ss)_size\\.name(.*)"
  ) %>%
  filter(q < p_threshold & passed_ss == TRUE) %>%   
  distinct(OTU) %>%
  pull(OTU)

# Define matrix of log-fold change data per size
data_mat <- output$res %>% 
  rename(OTU = taxon) %>%
  left_join(., taxonomy, by = "OTU") %>%
  # Keep DA_taxa
  filter(OTU %in% DA_taxa) %>%
  # change value to zero if not significant
  mutate(
    M   = ifelse(q_size.nameM   < p_threshold & passed_ss_size.nameM,   lfc_size.nameM, 0),
    L   = ifelse(q_size.nameL   < p_threshold & passed_ss_size.nameL,   lfc_size.nameL, 0),
    XL  = ifelse(q_size.nameXL  < p_threshold & passed_ss_size.nameXL,  lfc_size.nameXL, 0),
    XXL = ifelse(q_size.nameXXL < p_threshold & passed_ss_size.nameXXL, lfc_size.nameXXL, 0)
  ) %>%
  dplyr::select(Genus, M, L, XL, XXL) %>% # Species - new name
  rename_with(~ paste0(., "-S"), c("M", "L", "XL", "XXL")) %>%
  column_to_rownames(var = "Genus") %>%
  as.matrix()

# Metabolism
m_df <- as.data.frame(data_mat) %>%
  rownames_to_column(var = "Genus") %>%
  get_metabolism(., metab_fname) 

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

ht <- Heatmap(
  data_mat,
  # columns
  column_title = NULL, #"Differential Abundance (log fold change)", 
  cluster_columns = FALSE, # changes sample order
  column_names_rot = 0,
  column_names_centered = TRUE,
  # heatmap legend
  heatmap_legend_param = list(
    title = NULL, #"log fold change", 
    direction = "horizontal"
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
    height = 8, # height in inches; can adjust
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "top") #, annotation_legend_side = "top") 
draw(lgd, x = unit(0.66, "npc"), y = unit(0.98, "npc"), just = c("right", "top"))
dev.off()

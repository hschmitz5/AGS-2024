rm(list = ls())
source("./code/R/01_load_ps.R")
source("./code/R/02_metab_and_DA.R")
library(writexl)
library(ComplexHeatmap)
library(circlize) # for colorRamp2

write2excel <- FALSE
n_display_low <- 30

fname_excel <- "./data/ANCOM_ASV.xlsx"

fname_high  <- "./figures/DA_ancom_high.png"
fname_low   <- "./figures/DA_ancom_low.png"

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.6 # same as cell_h

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

rel_ab_cutoff <- 0.5

# ------ Process Data -----

# load differential abundance data
output <- readRDS(ancom_fname)

# define taxa in which at least one sample has abundance > rel_ab_cutoff
high_ab_taxa <- get_rel(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(OTU) %>%
  pull(OTU)

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
process_lfc <- function(df, taxa) {
  df %>%
    filter(OTU %in% taxa) %>%
    # change value to zero if not significant
    mutate(
      M   = ifelse(q_size.nameM   < p_threshold & passed_ss_size.nameM,   lfc_size.nameM, 0),
      L   = ifelse(q_size.nameL   < p_threshold & passed_ss_size.nameL,   lfc_size.nameL, 0),
      XL  = ifelse(q_size.nameXL  < p_threshold & passed_ss_size.nameXL,  lfc_size.nameXL, 0),
      XXL = ifelse(q_size.nameXXL < p_threshold & passed_ss_size.nameXXL, lfc_size.nameXXL, 0)
    ) %>%
    dplyr::select(Genus, M, L, XL, XXL) %>% # Species - new name
    mutate(
      mean_lfc = rowMeans(across(M:XXL, abs), na.rm = TRUE)
    ) %>%
    rename_with(~ paste0(., "-S"), c("M", "L", "XL", "XXL")) %>%
    arrange(desc(mean_lfc))
}

taxonomy <- get_taxonomy(ps)

# Read in data and add new names
res_prim <- output$res %>% 
  rename(OTU = taxon) %>%
  left_join(., taxonomy, by = "OTU")

# Apply function to high and low abundance taxa
lfc_high <- process_lfc(res_prim, DA_taxa) 
  
# lfc_low  <- process_lfc(res_prim, DA_taxa$low_ab)

### For plotting
fig_high <- lfc_high %>%
  dplyr::select(-mean_lfc) %>%
  tibble::column_to_rownames("Genus") %>%
  as.matrix()

# ---- Plotting

create_heatmap <- function(mat, rowname_w = NULL, col_title = NULL) {
  n_cols <- ncol(mat)
  n_rows <- nrow(mat)
  
  # Labels
  row_labels <- rownames(mat)
  # italicize species + _g, but not _f/_o/_c/_p
  italic_rows <- !grepl("_(f|o|c|p)(?:_|$|-)", row_labels)
  row_fontface <- ifelse(italic_rows, "italic", "plain")
  
  args <- list(
    mat,
    name = "log fold\nchange",
    cluster_columns = FALSE,
    show_heatmap_legend = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = col_title,
    #col = ht_colors,
    # size
    width  = unit(n_cols * cell_w, "inches"),
    height = unit(n_rows * cell_h, "inches"),
    row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
    column_names_gp = gpar(fontsize = col_fontsize)
  )
  # Only add row_names_max_width if it is not NULL
  if (!is.null(rowname_w)) {
    args$row_names_max_width <- unit(rowname_w, "inches")
  }
  
  do.call(Heatmap, args)
}




ht_high <- create_heatmap(fig_high, NULL, paste0("ANCOM-BC2 DA Taxa"))

col_fun <- colorRamp2(c(min(fig_high), 0, max(fig_high)), c("blue", "white", "red"))

lgd <- Legend(
  col_fun = col_fun,
  title = "log fold change",
  direction = "horizontal",
  title_position = "topcenter",
  at = pretty(c(min(fig_high), max(fig_high))),
  labels = pretty(c(min(fig_high), max(fig_high)))
)

# Define rowname width in first plot
ht_grob <- draw(ht_high, show_heatmap_legend = FALSE)
ht1 <- ht_grob@ht_list[[1]]
fig_props <- ht1@layout$layout_size

rowname_width_in <- convertWidth(
  fig_props$row_names_right_width,
  "inches", valueOnly = TRUE
)

# Save image
common_width <- ncol(fig_high) * cell_w + rowname_width_in + 2 # in
png(fname_high, 
    width = common_width,
    height = nrow(fig_high) * cell_h + 2,
    units = "in", res = 300)
draw(ht_high)
draw(lgd, x = unit(0.41, "npc"), y = unit(0.99, "npc"), just = c("center", "top"))
dev.off()
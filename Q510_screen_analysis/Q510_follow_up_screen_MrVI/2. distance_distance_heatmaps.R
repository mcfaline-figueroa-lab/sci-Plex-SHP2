library(monocle3)
library(tidyverse)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# ================================================================================
# distance-distance heatmap remade in R
# ================================================================================

heatmap_mat <- read_csv("SHP2_Q510_24hr_NoEGF_mrvi_dist_df_attention.csv", 
                        col_names = T) %>%
  dplyr::rename(X1 = 1) %>%
  column_to_rownames(var = "X1")

hmcols <- viridis::plasma(n = 50)
hm <- ComplexHeatmap::Heatmap(matrix = heatmap_mat %>% as.matrix(),
                              cluster_columns = T, cluster_rows = T, 
                              clustering_method_columns = "ward.D2", 
                              clustering_method_rows = "ward.D2",
                              col = hmcols,
                              name = "Distance", 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          legend_height = unit(0.1, "cm"),
                                                          legend_width = unit(2, "cm"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          # at = c(0, 5, 10),
                                                          # labels = c("0", "1", "2", "3", "4", "5"),
                                                          # labels = c("<0", "5", ">10"), 
                                                          direction = "horizontal"),
                              column_split = 3,
                              row_split = 3,
                              column_gap = unit(0.1, "line"),
                              row_gap = unit(0.1, "line"),
                              column_title = NULL,
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 6),
                              column_names_gp = gpar(fontsize = 5),
                              row_title = NULL,
                              show_column_names = F,
                              show_row_names = T,
                              row_names_gp = gpar(fontsize = 4), 
                              row_names_side = "right",
                              row_dend_gp = gpar(lwd = 0.5),
                              column_dend_gp = gpar(lwd = 0.5),
                              show_heatmap_legend = T)

ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "right",
                     legend_grouping = "original")


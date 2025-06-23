library(tidyverse)
library(monocle3)
library(kableExtra)
library(ComplexHeatmap)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

treatment_diff_test <- readRDS("SHP2_mut_EGF_treatment_diff_test.RDS")
treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

DEG_space <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(q_value <= 0.05) %>%
  distinct(id) %>%
  mutate(signif = "Yes")

cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

colData(cds)$EGF_dose_group <- case_when(
  colData(cds)$EGF_dose == "0" ~ "No EGF",
  colData(cds)$EGF_dose %in% c("12.5", "25", "50") ~ "Low",
  TRUE ~ "High"
)
colData(cds)$EGF_dose_group = factor(colData(cds)$EGF_dose_group, 
                                     levels = c("No EGF", "Low", "High"))

cds_subset <- cds[DEG_space$id,
                  colData(cds)$timepoint == "24hr" &
                    colData(cds)$SHP2_mut %in% c("KO", "WT")]

# =============================================================================
# Gene Modules
# =============================================================================

gene_module_df <- find_gene_modules(cds_subset[DEG_space$id,], resolution=1e-2) %>%
  left_join(rowData(cds_subset) %>% as.data.frame())

ggplot(gene_module_df, aes(x = dim_1, y = dim_2)) +
  geom_point(aes(color = module)) +
  monocle3:::monocle_theme_opts()
ggsave("UMAP_gene_module_DEG_EGF_WT_vs_KO.png", dpi = 600, height = 3, width = 3.5)

colData(cds_subset)$EGF_dose_group <- case_when(
  colData(cds_subset)$EGF_dose == "0" ~ "No EGF",
  colData(cds_subset)$EGF_dose %in% c("12.5", "25", "50") ~ "Low",
  colData(cds_subset)$EGF_dose %in% c("100", "250", "500", "1000") ~ "High"
)
colData(cds_subset)$EGF_dose_group = factor(colData(cds_subset)$EGF_dose_group, 
                                     levels = c("No EGF", "Low", "High"))

cell_group_df <- colData(cds_subset) %>% as.data.frame() %>%
  select(cell = cell_ID,
         SHP2_mut,
         EGF_dose_group) %>%
  unite("Group", SHP2_mut, EGF_dose_group, sep = "_")
agg_mat <- aggregate_gene_expression(cds = cds_subset, 
                                     gene_group_df = gene_module_df,  
                                     cell_group_df = cell_group_df , 
                                     cell_agg_fun = "mean", 
                                     scale_agg_values = F, 
                                     gene_agg_fun = "sum") %>% as.matrix()
agg_mat <- t(scale(t(agg_mat)))
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))

col_anno_df <- data.frame("temp" = colnames(agg_mat)) %>%
  separate(col = temp, into = c("SHP2", "EGF Dose"), sep = "_", remove = F) %>%
  tibble::column_to_rownames("temp") %>%
  mutate(`EGF Dose` = factor(`EGF Dose`, levels = c("No EGF", "Low", "High")))

col_anno <- HeatmapAnnotation(df = col_anno_df,
                              col = list(`SHP2` = c("KO" = "#F8766D",
                                                    "WT" = "#00BFC4"),
                                         `EGF Dose` = c("No EGF" =  viridis::magma(n = 4)[2], # "grey90",
                                                        "Low" = viridis::magma(n = 4)[3],
                                                        "High" = viridis::magma(n = 4)[4])), 
                              show_legend = T,
                              annotation_legend_param = list(`SHP2` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                                           labels_gp = gpar(fontsize = 6)),
                                                             `EGF Dose` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                                               labels_gp = gpar(fontsize = 6))), 
                              annotation_name_gp = gpar(fontsize = 6)
)

row_anno_df <- gene_module_df %>%
  group_by(module) %>%
  dplyr::summarise(`# Genes` = n()) %>%
  mutate(module = paste0("Module ",module)) %>%
  tibble::column_to_rownames("module")

row_anno_df <- row_anno_df[rownames(agg_mat), ]

row_bar <- rowAnnotation(`# Genes` = anno_barplot(row_anno_df, 
                                                  bar_width = 0.4,
                                                  add_numbers = TRUE,
                                                  numbers_gp = gpar(fontsize = 6),
                                                  numbers_rot = 0,
                                                  numbers_offset = unit(0.5,'mm'),
                                                  width = unit(1.5,'cm'),
                                                   axis_param = list(gp = gpar(fontsize = 6),
                                                                     side = "bottom",
                                                                     labels_rot = 45,
                                                                     at = c(0,100,200,300),
                                                                     labels = c("0","100", "200", "300"))
), 
annotation_legend_param = list(fontsize = 6),
annotation_name_gp = gpar(fontsize = 6),
annotation_name_side = "bottom")

hmcol <- circlize::colorRamp2(c(2,1,0,-1,-2), RColorBrewer::brewer.pal(n = 5, "RdYlBu"))
hm <- ComplexHeatmap::Heatmap(matrix = agg_mat,
                              name = "Z-Scored\nMean Epxression",
                              # col = hmcol,
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                          labels_gp = gpar(fontsize = 6),
                                                          border = TRUE),
                              clustering_method_columns = "ward.D2", 
                              clustering_method_rows = "ward.D2",
                              top_annotation = col_anno,
                              right_annotation = row_bar,
                              column_title = "SHP2_EGF Group",
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                              row_title = "Gene Module",
                              row_title_side = "left",
                              row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                              show_column_names = F,
                              show_row_names = T,
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                              cluster_rows = T,
                              cluster_columns = T,
                              row_split = 2,
                              column_split = 3,
                              row_dend_width = unit(0.2,"cm"),
                              column_dend_height = unit(0.4,"cm"),
                              use_raster = FALSE
)

draw(hm)

pdf("gene_module_EGF_DEGs_WT_KO_heatmap.pdf", width = 5, height = 2, compress = F,
    bg = "white")
ht <- draw(hm, merge_legends = F); gene_sets <- gene_sets <- row_order(ht); col_clusters <- column_order(ht)
dev.off()

# write_csv(gene_module_df, "gene_module_results_EGF_DEGs_20250114.csv")


# ================================================================================
# Further characterize gene modules
# ================================================================================
gene_module_df <- read_csv("gene_module_results_EGF_DEGs_20250114.csv")

# Making silhouette score plots for gene module membership
dist.matrix <- dist(x = gene_module_df[,4:5])
sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = gene_module_df$module)), dist = dist.matrix)
gene_module_df$sil <- sil[,3]
mean_silhouette_score <- mean(gene_module_df$sil)

gene_module_df %>%
  mutate(barcode = rownames(.)) %>%
  arrange(module,-sil) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  mutate(module = paste0("Module ",module)) %>%
  ggplot() +
    geom_col(aes(barcode, sil, fill = module), width = 1, show.legend = FALSE) +
    geom_hline(yintercept = mean_silhouette_score, color = 'red', linetype = 'dashed') +
    annotate(geom = "text", x = 100, y = -0.25, size = 2,
             label = paste0("mean = ",round(mean_silhouette_score, 3))) +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette score') +
    viridis::scale_fill_viridis(discrete = TRUE) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
ggsave(paste0("silhouette_scores_gene_module_results_EGF_DEGs_20250114.png"), 
       dpi = 600, height = 2, width = 2.75)


# Hypergeometric GSA tests on modules
library(piano)
source("loadGSCSafe.R")

pathways.hallmark <- loadGSCSafe("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/h.all.v6.0.symbols.gmt")
pathways.reactome <- loadGSCSafe("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c2.cp.reactome.v2024.1.Hs.symbols.gmt")

pathways.all <- list("gsc" = c(pathways.hallmark[["gsc"]], pathways.reactome[["gsc"]]),
                     "addInfo" = c(pathways.hallmark[["addInfo"]], pathways.reactome[["addInfo"]]))

class(pathways.all) <- "GSC"

gene_universe <- treatment_diff_test %>%
  distinct(gene_short_name) %>%
  pull()

GSA_full <- data.frame()
for (mod in 1:4) {
  top_genes <- gene_module_df %>%
    filter(module == mod) %>%
    select(gene_short_name) %>%
    pull()
  
  test <- piano::runGSAhyper(genes = top_genes,
                             gsc = pathways.hallmark,
                             adjMethod = "none",
                             universe = gene_universe)
  GSAhyper_df <- as.data.frame(test$p.adj) %>%
    rownames_to_column(var = "gene_set")
  colnames(GSAhyper_df) <- c("gene_set","p_value")
  
  GSAhyper_df <- GSAhyper_df %>%
    mutate(module = mod)
  
  GSA_full <- bind_rows(GSA_full, GSAhyper_df)
  GSA_full$q_value <- p.adjust(GSA_full$p_value, method = "BH")
  
}


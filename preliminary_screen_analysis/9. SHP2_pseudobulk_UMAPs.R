library(tidyverse)
library(monocle3)

this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

colData(cds)$SHP2_mut <- factor(colData(cds)$SHP2_mut, levels = c("WT", "KO",
                                                                  "T42A", "T52S",
                                                                  "E76K", "R138Q",
                                                                  "E139D", "Y279C",
                                                                  "T468M", "T507K",
                                                                  "Q510E", "Q510K"))

# load SHP2 DEGs for both 24hr and 96hr timepoints
DEG_union <- readRDS("SHP2_time_treatment_diff_test.RDS")
DEG_union$q_value <- p.adjust(DEG_union$p_value, method = "BH")
DEG_union_SHP2 <- DEG_union %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(q_value <= 0.05) %>%
  distinct(id) %>%
  pull()

# =====================================================================================
# Pseudo-bulking by timepoint, mutant, EGF dose
# =====================================================================================

dir.create("pseudobulked")

# Aggregate counts by drug, dose and batch by summing the gene expression counts from nuclei 
# in each condition
aggregated_counts_by_drug_dose =
  colData(cds) %>%
  as.data.frame() %>% 
  mutate(condition = paste0(timepoint,"_",SHP2_mut,"_",EGF_dose)) %>%
  group_by(condition) %>%
  nest() %>%
  mutate(aggregated_counts = 
           purrr::map(.x = data,
                      .f = function(coldata_subset,
                                    cds){
                        counts(cds)[,coldata_subset$Cell] %>%
                          Matrix::rowSums()
                      },cds))

aggregated_count_matrix = 
  do.call(cbind, aggregated_counts_by_drug_dose$aggregated_counts)

colnames(aggregated_count_matrix) <- aggregated_counts_by_drug_dose$condition

pseudobulk_cds <- new_cell_data_set(expression_data = aggregated_count_matrix,
                                    cell_metadata = data.frame(row.names = colnames(aggregated_count_matrix), cell = colnames(aggregated_count_matrix)),
                                    gene_metadata = data.frame(row.names = row.names(aggregated_count_matrix), id = row.names(aggregated_count_matrix)))


rowData(pseudobulk_cds)$gene_short_name <- rowData(cds)[row.names(pseudobulk_cds),]$gene_short_name

colData(pseudobulk_cds)$timepoint <- sapply(colData(pseudobulk_cds)$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][1]})
colData(pseudobulk_cds)$SHP2_mut <- sapply(colData(pseudobulk_cds)$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][2]})
colData(pseudobulk_cds)$EGF_dose <- sapply(colData(pseudobulk_cds)$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][3]})

pseudobulk_cds <- detect_genes(pseudobulk_cds, min_expr = 0)

normalized_count_matrix <- normalized_counts(pseudobulk_cds[rowSums(exprs(pseudobulk_cds))>0,], 
                                             norm_method = "size_only")
normalized_count_matrix_long <- reshape2::melt(normalized_count_matrix %>% as.matrix)
colnames(normalized_count_matrix_long) <- c("id","cell","expression")

new_metadata <- data.frame(cell = unique(normalized_count_matrix_long$cell),
                           condition = unique(normalized_count_matrix_long$cell))

new_metadata$timepoint <- sapply(new_metadata$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][1]})
new_metadata$SHP2_mut <- sapply(new_metadata$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][2]})
new_metadata$EGF_dose <- sapply(new_metadata$cell,function(x){stringr::str_split(x, pattern = "_")[[1]][3]})

normalized_count_matrix_long <- left_join(normalized_count_matrix_long, new_metadata, by = "cell")

normalized_count_matrix_long <- normalized_count_matrix_long %>%
  group_by(id, timepoint) %>%
  dplyr::mutate(mean_log2_fc = log2((expression + 1)/(expression[SHP2_mut == "KO" & EGF_dose == 0] + 1)))

normalized_l2fc_matrix <- normalized_count_matrix_long %>%
  ungroup() %>%
  filter(SHP2_mut != "KO") %>%
  filter(id %in% DEG_union_SHP2) %>%
  mutate(cell = paste0(timepoint,"_",SHP2_mut,"_",EGF_dose)) %>%
  dplyr::select(id, cell, mean_log2_fc) %>%
  tidyr::spread(key = cell, value = mean_log2_fc) %>%
  as.data.frame()

row.names(normalized_l2fc_matrix) <- normalized_l2fc_matrix$id
normalized_l2fc_matrix$id <- NULL

temp_gene_meta <- data.frame(row.names = row.names(normalized_l2fc_matrix), 
                             id = row.names(normalized_l2fc_matrix)) %>%
  left_join(rowData(cds) %>% as.data.frame(),
            by = c("id")) %>%
  mutate(temp_id = id) %>%
  tibble::column_to_rownames(var = "temp_id")

temp_cell_meta <- data.frame(row.names = colnames(normalized_l2fc_matrix), 
                             cell = colnames(normalized_l2fc_matrix)) %>%
  separate(cell, into = c("timepoint", "SHP2_mut", "EGF_dose"), sep = "_", remove = F)

pseudobulk_log2fc_cds <- new_cell_data_set(expression_data = normalized_l2fc_matrix %>% as.matrix(),
                                           cell_metadata = temp_cell_meta,
                                           gene_metadata = temp_gene_meta)
colData(pseudobulk_log2fc_cds)$Size_Factor <- 1


pseudobulk_log2fc_cds <- preprocess_cds(pseudobulk_log2fc_cds, 
                                   method = "PCA", 
                                   num_dim = 5, 
                                   norm_method = "none", 
                                   scaling = F)

plot_pc_variance_explained(pseudobulk_log2fc_cds)

# pseudobulk_log2fc_cds <- align_cds(pseudobulk_log2fc_cds, 
#                                    method = "PCA",
#                                    alignment_group = "timepoint", 
#                                    alignment_k = 10,
#                                    verbose = T)

pseudobulk_log2fc_cds <- reduce_dimension(pseudobulk_log2fc_cds, 
                                     reduction_method = "UMAP", 
                                     preprocess_method = "PCA",
                                     umap.n_neighbors = 10, 
                                     umap.min_dist = 0.2, 
                                     max_components = 2,
                                     verbose = TRUE)

colData(pseudobulk_log2fc_cds)$test_col <- case_when(
  colData(pseudobulk_log2fc_cds)$SHP2_mut == "R138Q" ~ "R138Q",
  TRUE ~ "Other"
)

colData(pseudobulk_log2fc_cds)$test_col <- factor(
  colData(pseudobulk_log2fc_cds)$test_col,
  levels = c("R138Q", "Other")
)

colData(pseudobulk_log2fc_cds)$EGF_dose <- factor(
  colData(pseudobulk_log2fc_cds)$EGF_dose,
  levels = c("0", "12.5", "25", "50", "100", "250", "500", "1000")
)

colData(pseudobulk_log2fc_cds)$test_col <- case_when(
  colData(pseudobulk_log2fc_cds)$EGF_dose == "0" ~ "0ng/ul",
  TRUE ~ "Other"
)

colData(pseudobulk_log2fc_cds)$test_col <- factor(
  colData(pseudobulk_log2fc_cds)$test_col,
  levels = c("0ng/ul", "Other")
)

plot_cells(pseudobulk_log2fc_cds, color_cells_by = "timepoint", 
           cell_size = 1.5, 
           alpha = 1, 
           cell_stroke = 0.01,
           label_cell_groups = F) + 
  # facet_wrap(~SHP2_mut) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  viridis::scale_color_viridis(name = "Timepoint", discrete = TRUE) +
  guides(color=guide_legend(title="Timepoint", byrow = TRUE, override.aes = list(size=2, alpha = 1)))
ggsave("pseudobulked/UMAP_timepoint_without_align.png",
       dpi = 600, height = 1.5, width = 2)

colData(pseudobulk_log2fc_cds)$EGF_group <- case_when(
  colData(pseudobulk_log2fc_cds)$EGF_dose == "0" ~ "No EGF",
  colData(pseudobulk_log2fc_cds)$EGF_dose %in% c("12.5", "25", "50") ~ "Low",
  TRUE ~ "High"
)

plot_cells(pseudobulk_log2fc_cds, color_cells_by = "EGF_group", 
           cell_size = 1.5, 
           alpha = 1, 
           cell_stroke = 0.01,
           label_cell_groups = F) + 
  # facet_wrap(~SHP2_mut) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  viridis::scale_color_viridis(name = "EGF Group", discrete = TRUE, option = "magma", direction = -1) +
  guides(color=guide_legend(title="EGF Group", byrow = TRUE, override.aes = list(size=2, alpha = 1)))
ggsave("pseudobulked/UMAP_group_without_align.png",
       dpi = 600, height = 1.5, width = 2)

plot_cells(pseudobulk_log2fc_cds, color_cells_by = "test_col", 
           cell_size = 1.5, 
           alpha = 1, 
           cell_stroke = 0.01,
           label_cell_groups = F) + 
  # facet_wrap(~SHP2_mut) +
  theme(legend.position = "right",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  viridis::scale_color_viridis(name = "R138Q", discrete = TRUE, option = "magma", direction = -1) +
  guides(color=guide_legend(title="R138Q", byrow = TRUE, override.aes = list(size=2, alpha = 1)))
ggsave("pseudobulked/UMAP_R138Q_without_align.png",
       dpi = 600, height = 1.5, width = 2)

saveRDS(pseudobulk_log2fc_cds,
        file = "pseudobulked/SHP2_pseudobulk_log2FC_to_KO-EGF0_cds.RDS")

pseudobulk_log2fc_cds <- readRDS("pseudobulked/SHP2_pseudobulk_log2FC_to_KO-EGF0_cds.RDS")

colData(pseudobulk_log2fc_cds)$UMAP1 <- reducedDims(pseudobulk_log2fc_cds)[["UMAP"]][,1]
colData(pseudobulk_log2fc_cds)$UMAP2 <- reducedDims(pseudobulk_log2fc_cds)[["UMAP"]][,2]

col_data_for_printing <- colData(pseudobulk_log2fc_cds) %>% as.data.frame() %>%
  select(-Size_Factor, -test_col)

write_csv(col_data_for_printing,
          col_names = T,
          file = "pseudobulked/SHP2_pseudobulk_log2fc_cds_UMAP_coordinates_no_time_align.csv")


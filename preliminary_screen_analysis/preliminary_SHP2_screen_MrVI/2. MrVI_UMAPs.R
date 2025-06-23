library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)
source("calculate_aggreg_expression.R")

# First -- run chunk at bottom of sheet to establish modified plot_cells function

cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.05,]$id

cds <- cds[,colData(cds)$timepoint == "24hr"]

temp_col_data_rownames <- rownames(colData(cds) %>% as.data.frame())

# ================================================================================
# UMAPs of U latent spaces
# ================================================================================

MrVI_u_space <- read_csv("SHP2_24hr_MrVI_u_factor.csv",
                         col_names = F) %>%
  as.data.frame()

row.names(MrVI_u_space) <- temp_col_data_rownames

reducedDims(cds)[["PCA"]] <- MrVI_u_space

cds <- reduce_dimension(cds, 
                        preprocess_method = "PCA", 
                        verbose = T, umap.n_neighbors = 25,
                        umap.min_dist = 0.2
                        )

# reducedDims(cds)[["UMAP"]]

plot_cells(cds,
           color_cells_by = "SHP2_mut",
           cell_size = .5, alpha = 0.8, 
           label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # facet_grid(hash_plate~SHP2_mut) +
  # viridis::scale_color_viridis(option = "magma", discrete = T)
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))
ggsave("Uspace_UMAP_SHP2.png",
       dpi = 900, height = 3, width = 3)

colData(cds)$EGF_dose <- factor(colData(cds)$EGF_dose,
                                levels = c("0", "12.5", "25", "50", "100", "250", "500", "1000"))

plot_cells(cds,
           color_cells_by = "EGF_dose",
           cell_size = .5, alpha = 0.8, 
           label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # facet_grid(hash_plate~SHP2_mut) +
  viridis::scale_color_viridis("EGF (ng/mL)",option = "magma", discrete = T)
# scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))
ggsave("Uspace_UMAP_EGFdose.png",
       dpi = 900, height = 3, width = 3)

plot_cells(cds, reduction_method = "UMAP",
           color_cells_by = "proliferation_index",
           cell_size = .5, alpha = 1, 
           label_cell_groups = F) +
  # facet_wrap(~log10_dose) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # guides(color = "none") +
  viridis::scale_color_viridis(name = "Proliferation\nIndex", discrete = FALSE, option = "D")
ggsave("Uspace_UMAP_proliferation_index.png",
       dpi = 900, height = 3, width = 3)

plot_cells(cds, reduction_method = "UMAP",
           genes = "PTPN11",
           scale_to_range = F,
           cell_size = .5, alpha = 1, 
           label_cell_groups = F) +
  # facet_wrap(~log10_dose) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # guides(color = "none") +
  viridis::scale_color_viridis(name = "PTPN11\nExpression", discrete = FALSE, option = "D")
ggsave("Uspace_UMAP_SHP2expression.png",
       dpi = 900, height = 3, width = 3)

colData(cds)$UMAP_1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP_2 <- reducedDims(cds)[["UMAP"]][,2]

colData(cds)$PTPN11_epxression <- calculate_aggregate_expression_score(cds, signature_genes = "PTPN11")

col_data_for_printing <- colData(cds) %>% as.data.frame() %>%
  select(UMAP_1, UMAP_2, timepoint, SHP2_mut, EGF_dose, PTPN11_epxression, proliferation_index)

write_csv(x = col_data_for_printing, 
          file = "u-space_allcells.csv",
          col_names = T)

# ================================================================================
# UMAPs of Z latent spaces
# ================================================================================

MrVI_z_space <- read_csv("SHP2_24hr_MrVI_z_factor.csv",
                         col_names = F) %>%
  as.data.frame()

row.names(MrVI_z_space) <- temp_col_data_rownames

reducedDims(cds)[["PCA"]] <- MrVI_z_space

cds <- reduce_dimension(cds, 
                        preprocess_method = "PCA", 
                        umap.min_dist = 0.05, umap.n_neighbors = 15,
                        verbose = T)

plot_cells(cds,
           color_cells_by = "SHP2_mut",
             cell_size = .5, alpha = 0.8, 
             label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # facet_grid(hash_plate~SHP2_mut) +
  # viridis::scale_color_viridis(option = "magma", discrete = T)
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))
ggsave("Zspace_UMAP_SHP2.png",
       dpi = 900, height = 3, width = 3)

colData(cds)$EGF_dose <- factor(colData(cds)$EGF_dose, levels = c("0", "12.5", "25", "50", "100", "250", "500", "1000"))

plot_cells(cds,
           color_cells_by = "EGF_dose",
           cell_size = .5, alpha = 0.8, 
           label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # facet_grid(hash_plate~SHP2_mut) +
  viridis::scale_color_viridis("EGF (ng/mL)",option = "magma", discrete = T)
  # scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))
ggsave("Zspace_UMAP_EGFdose.png",
       dpi = 900, height = 3, width = 3)

plot_cells(cds, reduction_method = "UMAP",
             genes = "PTPN11",
           scale_to_range = F,
             cell_size = .5, alpha = 1, 
             label_cell_groups = F) +
  # facet_wrap(~log10_dose) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # guides(color = "none") +
  viridis::scale_color_viridis(name = "Expression", discrete = FALSE, option = "D")
ggsave("Zspace_UMAP_SHP2expression.png",
       dpi = 900, height = 3, width = 3)

colData(cds)$plot_color <- case_when(
  colData(cds)$SHP2_mut == "R138Q" ~ "R138Q",
  colData(cds)$EGF_dose == "0" & !colData(cds)$SHP2_mut %in% c("Q510E", "KO")  ~ "No EGF",
  TRUE ~ NA
)
plot_cells(cds,
           color_cells_by = "plot_color",
           # alpha_cells_by = "plot_color",
           label_cell_groups = F,
           cell_stroke = 0.1, cell_size = 1
           ) +
  scale_alpha_manual(values = c(1,1), na.value = 0.5) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))

# Clustering in the z-space

cds <- cluster_cells(cds,
                     resolution = 1.5e-3,
                     reduction_method = "PCA",
                     num_iter = 1,
                     random_seed = 1995L)

colData(cds)$Cluster <- clusters(cds, reduction_method = "PCA")
length(unique(colData(cds)$Cluster))
colData(cds)$Partition <- partitions(cds, reduction_method = "PCA")
length(unique(colData(cds)$Partition))

plot_cells(cds, color_cells_by = "Cluster", cell_size = .5, alpha = 0.8, label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 11, "Set1")) +
  guides(color=guide_legend(title="Cluster", byrow = TRUE, nrow = 1,override.aes = list(size=2)))
ggsave("Zspace_UMAP_cluster.png",
       dpi = 600, height = 3, width = 3)

test_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(SHP2_mut, EGF_dose, Cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(SHP2_mut, EGF_dose) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(data = test_plot, aes(x = EGF_dose)) +
  facet_wrap(~SHP2_mut) +
  geom_bar(stat = "identity", aes(fill= Cluster, y = prop), position = "stack") +
  # geom_bar(aes(fill = Cluster, y = after_stat(prop)), width = 0.75) +
  scale_y_continuous(labels = scales::percent) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, "Set1")) +
  guides(fill=guide_legend(title="Cluster", override.aes = list(size=2), ncol = 1)) +
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.margin = unit(c(0,0,0,0), "cm"),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percentage") +
  xlab("EGF Dose (ng/mL)")
ggsave("Zspace_cluster_breakdown_SHP2_mut.png",
       dpi = 600, height = 4, width = 4.5)

cluster_percent_mat <- colData(cds) %>% as.data.frame() %>%
  group_by(SHP2_mut, EGF_dose, Cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(SHP2_mut, EGF_dose) %>%
  mutate(prop = count/sum(count)) %>%
  select(SHP2_mut, EGF_dose, Cluster, prop) %>%
  unite(col = "id", sep = "_", SHP2_mut, EGF_dose) %>%
  pivot_wider(id_cols = id, names_from = Cluster, values_from = prop, values_fill = 0) %>%
  column_to_rownames("id") %>%
  t()

cluster_percent_mat <- cluster_percent_mat[c("1", "2", "3", "4", "5", "6"),]

gene_meta <- data.frame(id = c("1", "2", "3", "4", "5", "6"), 
                        gene_short_name = NA)
row.names(gene_meta) <- gene_meta$id
cell_meta <- data.frame(id = colnames(cluster_percent_mat), id2 = colnames(cluster_percent_mat)) %>%
  separate(id2, into = c("SHP2_mut", "EGF_dose"), sep = "_")
row.names(cell_meta) <- cell_meta$id


cds_corr <- new_cell_data_set(cluster_percent_mat %>% as.matrix(),
                              gene_metadata = gene_meta,
                              cell_metadata = cell_meta)

reducedDims(cds_corr)[["PCA"]] <- t(exprs(cds_corr))

#plot_pc_variance_explained(cds)

cds_corr <- reduce_dimension(cds_corr, umap.min_dist = 0.25, umap.n_neighbors = 15)
colData(cds_corr)$UMAP1 <- reducedDims(cds_corr)[["UMAP"]][,1]
colData(cds_corr)$UMAP2 <- reducedDims(cds_corr)[["UMAP"]][,2]

plot_cells(cds_corr, color_cells_by = "SHP2_mut", cell_size = 0, alpha = 0.8, label_cell_groups = F) +
  geom_point(aes(color = SHP2_mut, size = EGF_dose)) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, "Set3"))
  # guides(color=guide_legend(title="Cluster", byrow = TRUE, nrow = 1,override.aes = list(size=2)))

library(ComplexHeatmap)

hmcol <- viridis::viridis(n = 20)
hm <- ComplexHeatmap::Heatmap(matrix = cluster_percent_mat, 
                              name = "mat",
                              col = hmcol,
                              heatmap_legend_param = list(title = "% of Cells\nin Cluster",
                                                          title_gp = gpar(fontface = "bold", fontsize = 6),
                                                          at = c(0, 0.5, 1),
                                                          labels = c("0", "50", "100"),
                                                          labels_gp = gpar(fontsize = 5),
                                                          grid_width = unit(0.3, "cm"),
                                                          legend_height = unit(1.2, "cm")),
                              # col = hmcols,
                              show_column_names = F, 
                              column_names_gp = gpar(fontsize = 6), column_names_rot = 45,
                              column_split = 8,
                              clustering_method_columns = "ward.D2",
                              # column_km = 5,
                              # row_split = 3,
                              column_title = NULL,
                              cluster_columns = TRUE, 
                              cluster_rows = TRUE,
                              row_title = "Cluster", 
                              row_title_gp = gpar(fontsize = 6),
                              show_row_names = T, 
                              row_names_gp = gpar(fontsize = 6), row_names_side = "right",
                              use_raster = F)

dir.create("z-space_hm")
png("z-space_hm/SHP2_mut_24hr_MrVI_z-space_cluster_hm.png",width=7,height=4,units="in",res=1200)
ht <- ComplexHeatmap::draw(hm, legend_grouping = "adjusted", merge_legends = T); col_clusters <- column_order(ht)
dev.off()

colData(cds)$UMAP_1 <- reducedDims(cds)[["UMAP"]][,1]
colData(cds)$UMAP_2 <- reducedDims(cds)[["UMAP"]][,2]

colData(cds)$PTPN11_epxression <- calculate_aggregate_expression_score(cds, signature_genes = "PTPN11")

col_data_for_printing <- colData(cds) %>% as.data.frame() %>%
  select(UMAP_1, UMAP_2, timepoint, SHP2_mut, EGF_dose, Cluster, PTPN11_epxression, proliferation_index)

write_csv(x = col_data_for_printing, 
          file = "z-space_allcells.csv",
          col_names = T)


# ================================================================================
# UMAPs of Z latent spaces without KO and transfection control
# ================================================================================
dir.create("exclude_KOcontrols")

MrVI_z_space <- read_csv("SHP2_24hr_MrVI_z_factor.csv",
                         col_names = F) %>%
  as.data.frame()

row.names(MrVI_z_space) <- temp_col_data_rownames

reducedDims(cds)[["PCA"]] <- MrVI_z_space

cds_subset <- reduce_dimension(cds[,!colData(cds)$SHP2_mut %in% c("KO", "Q510E")], 
                        preprocess_method = "PCA", 
                        umap.min_dist = 0.05, umap.n_neighbors = 10,
                        verbose = T)

plot_cells(cds_subset,
           color_cells_by = "SHP2_mut",
           cell_size = .5, alpha = 0.8, 
           label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 10, name = "Set3"))
ggsave("exlude_KOcontrols/Zspace_UMAP_SHP2.png",
       dpi = 900, height = 3, width = 3)

colData(cds_subset)$EGF_dose <- factor(colData(cds_subset)$EGF_dose, levels = c("0", "12.5", "25", "50", "100", "250", "500", "1000"))

plot_cells(cds_subset,
           color_cells_by = "EGF_dose",
           cell_size = .5, alpha = 0.8, 
           label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # facet_grid(hash_plate~SHP2_mut) +
  viridis::scale_color_viridis("EGF (ng/mL)",option = "magma", discrete = T)
# scale_color_manual(values = RColorBrewer::brewer.pal(n = 12, name = "Set3"))
ggsave("exclude_KOcontrols/Zspace_UMAP_EGFdose.png",
       dpi = 900, height = 3, width = 3)

plot_cells(cds_subset, reduction_method = "UMAP",
           genes = "PTPN11",
           scale_to_range = F,
           cell_size = .5, alpha = 1, 
           label_cell_groups = F) +
  # facet_wrap(~log10_dose) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = .01, r = .1, b = .1, l = .1, unit = "pt"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  # guides(color = "none") +
  viridis::scale_color_viridis(name = "Expression", discrete = FALSE, option = "D")
ggsave("exclude_KOcontrols/Zspace_UMAP_SHP2expression.png",
       dpi = 900, height = 3, width = 3)

cds_subset <- cluster_cells(cds_subset,
                            resolution = 1e-3,
                            reduction_method = "PCA",
                            num_iter = 5,
                            random_seed = 1995L)

colData(cds_subset)$Cluster <- clusters(cds_subset, reduction_method = "PCA")
length(unique(colData(cds_subset)$Cluster))

plot_cells(cds_subset, color_cells_by = "Cluster", cell_size = .5, alpha = 0.8, label_cell_groups = F) +
  theme(legend.position = "bottom",
        text = element_text(size = 6),
        legend.key.width = unit(0.3,"line"),
        legend.key.height = unit(0.3,"line"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt"),
        legend.spacing.y = unit(0.1, "cm"),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 11, "Set1")) +
  guides(color=guide_legend(title="Cluster", byrow = TRUE, nrow = 1,override.aes = list(size=2)))
ggsave("exclude_KOcontrols/Zspace_UMAP_cluster.png",
       dpi = 600, height = 3, width = 3)

test_plot <- colData(cds_subset) %>% as.data.frame() %>%
  dplyr::group_by(SHP2_mut, EGF_dose, Cluster) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(SHP2_mut, EGF_dose) %>%
  dplyr::mutate(prop = count / sum(count)) %>%
  ungroup()

ggplot(data = test_plot, aes(x = EGF_dose)) +
  facet_wrap(~SHP2_mut, nrow = 2) +
  geom_bar(stat = "identity", aes(fill= Cluster, y = prop), position = "stack") +
  # geom_bar(aes(fill = Cluster, y = after_stat(prop)), width = 0.75) +
  scale_y_continuous(labels = scales::percent) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, "Set1")) +
  guides(fill=guide_legend(title="Cluster", override.aes = list(size=2), ncol = 1)) +
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.margin = unit(c(0,0,0,0), "cm"),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percentage") +
  xlab("EGF Dose (ng/mL)")
ggsave("exclude_KOcontrols/Zspace_cluster_breakdown_SHP2_mut.png",
       dpi = 600, height = 2.5, width = 4.5)

test_plot <- colData(cds_subset) %>% as.data.frame() %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    # EGF_dose %in% c("50", "100", "250") ~ "Mid",
    TRUE ~ "High"
  )) %>%
  dplyr::group_by(SHP2_mut, EGF_dose_group, Cluster) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(SHP2_mut, EGF_dose_group) %>%
  dplyr::mutate(prop = count / sum(count)) %>%
  ungroup() %>%
  mutate(EGF_dose_group = factor(EGF_dose_group, levels = c("No EGF", "Low", "Mid", "High")))

test_plot <- test_plot %>%
  mutate(EGF_dose_group = factor(EGF_dose_group, levels = c("No EGF", "Low", "High")))

ggplot(data = test_plot, aes(x = EGF_dose_group)) +
  facet_wrap(~SHP2_mut, nrow = 2) +
  geom_bar(stat = "identity", aes(fill = as.factor(Cluster), y = prop), position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  monocle3:::monocle_theme_opts() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, "Set1")) +
  guides(fill=guide_legend(title="Cluster", override.aes = list(size=2), ncol = 1)) +
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.margin = unit(c(0,0,0,0), "cm"),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percentage") +
  xlab("EGF Dose (ng/mL)")
ggsave("exclude_KOcontrols/Zspace_cluster_breakdown_SHP2_mut_EGFbinned.png",
       dpi = 600, height = 2.5, width = 4.5)

colData(cds_subset)$UMAP_1 <- reducedDims(cds_subset)[["UMAP"]][,1]
colData(cds_subset)$UMAP_2 <- reducedDims(cds_subset)[["UMAP"]][,2]

colData(cds_subset)$PTPN11_epxression <- calculate_aggregate_expression_score(cds_subset, 
                                                                              signature_genes = "PTPN11")

colData(cds_subset)$EGF_dose_group <- case_when(
  colData(cds_subset)$EGF_dose == "0" ~ "No EGF",
  colData(cds_subset)$EGF_dose %in% c("12.5", "25", "50") ~ "Low",
  TRUE ~ "High"
)

col_data_for_printing <- colData(cds_subset) %>% as.data.frame() %>%
  select(UMAP_1, UMAP_2, timepoint, SHP2_mut, EGF_dose, EGF_dose_group, Cluster, PTPN11_epxression, proliferation_index)

write_csv(x = col_data_for_printing, 
          file = "exclude_KOcontrols/SHP2_MrVI_z-space_UMAPcoord_excludeKO.csv",
          col_names = T)

write_csv(x = test_plot,
          file = "exclude_KOcontrols/SHP2_MrVI_z-space_cluster_proportions_excludeKO.csv",
          col_names = T)

test_plot <- read_csv("Uexclude_KOcontrols/SHP2_MrVI_z-space_cluster_proportions_excludeKO.csv")

# saveRDS(cds_subset, file = "SHP2_cds_noKO_clustered_z-space.RDS")
# cds_subset <- readRDS("SHP2_cds_noKO_clustered_z-space.RDS")

col_data_sil <- col_data_for_printing
dist.matrix <- dist(x = col_data_sil[,1:2])
sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = col_data_sil$Cluster)), dist = dist.matrix)
col_data_sil$sil <- sil[,3]
mean_silhouette_score <- mean(col_data_sil$sil)

col_data_sil %>%
  mutate(barcode = rownames(.)) %>%
  arrange(Cluster,-sil) %>%
  mutate(barcode = factor(barcode, levels = barcode)) %>%
  ggplot() +
  geom_col(aes(barcode, sil, fill = Cluster), width = 1, show.legend = FALSE) +
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
ggsave(paste0("silhouette_scores_MrVI_z-space_noKO.png"), 
       dpi = 600, height = 2, width = 2.75)

test_plot <- colData(cds_subset) %>% as.data.frame() %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    TRUE ~ "High"
  )) %>%
  dplyr::group_by(SHP2_mut, EGF_dose_group, Cluster) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::mutate(proportion_of_cluster = count / sum(count)) %>%
  ungroup()

write_csv(x = test_plot,
          file = "UMAP_coordinates_for_Anne/exclude_KOcontrols/SHP2_MrVI_z-space_mutant_EGFdose_proportions_excludeKO.csv",
          col_names = T)

ggplot(data = test_plot, aes(x = Cluster)) +
  # facet_wrap(~SHP2_mut, nrow = 2) +
  geom_bar(stat = "identity", aes(fill= EGF_dose_group, y = proportion_of_cluster), position = "stack") +
  # geom_bar(aes(fill = Cluster, y = after_stat(prop)), width = 0.75) +
  scale_y_continuous(labels = scales::percent) +
  monocle3:::monocle_theme_opts() +
  # scale_fill_manual(values = RColorBrewer::brewer.pal(n = 8, "Set1")) +
  guides(fill=guide_legend(title="Cluster", override.aes = list(size=2), ncol = 1)) +
  theme(legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(0,0,0,0, "cm"),
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Percentage") +
  xlab("Cluster")

# ================================================================================
# DEG test of cluster 5 vs all other clusters
# ================================================================================

cds_subset <- readRDS("SHP2_cds_noKO_clustered_z-space.RDS")

colData(cds_subset)$ClusterDEG <- case_when(
  colData(cds_subset)$Cluster == 5 ~ "Cluster5",
  TRUE ~ "Other"
)

colData(cds_subset)$ClusterDEG <- factor(colData(cds_subset)$ClusterDEG, levels = c("Other", "Cluster5"))

treatment_diff_test <- fit_models(cds = cds_subset[expressed_genes,], 
                                  model_formula_str = "~ ClusterDEG",
                                  expression_family = "quasipoisson",
                                  cores = 1, 
                                  verbose = TRUE)
treatment_diff_test <- coefficient_table(treatment_diff_test) %>% dplyr::select(-model,-model_summary)
treatment_diff_test$SHP2_mut <- stringr::str_remove(treatment_diff_test$term,
                                                    "ClusterDEG")

treatment_diff_test$q_value <- p.adjust(p = treatment_diff_test$p_value, method = "BH")  

saveRDS(treatment_diff_test, "exclude_KOcontrols/Cluster5_full_treatment_diff_test.RDS")
treatment_diff_test <- readRDS("exclude_KOcontrols/Cluster5_full_treatment_diff_test.RDS")

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "ClusterDEG", term)) %>%
  select(Cluster = SHP2_mut, id, gene_short_name, normalized_effect, q_value)

readr::write_csv(file = "exclude_KOcontrols/Cluster5_vs_All_MrVI_DEGs.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)

# cluster 5 vs just Q510K and T507K that are distinctly absent
colData(cds_subset)$ClusterDEG <- case_when(
  colData(cds_subset)$Cluster == 5 ~ "Cluster5",
  colData(cds_subset)$SHP2_mut %in% c("Q510K", "T507K") ~ "Other",
  TRUE ~ NA
)

colData(cds_subset)$ClusterDEG <- factor(colData(cds_subset)$ClusterDEG, levels = c("Other", "Cluster5"))

treatment_diff_test <- fit_models(cds = cds_subset[expressed_genes,], 
                                  model_formula_str = "~ ClusterDEG",
                                  expression_family = "quasipoisson",
                                  cores = 1, 
                                  verbose = TRUE)
treatment_diff_test <- coefficient_table(treatment_diff_test) %>% dplyr::select(-model,-model_summary)
treatment_diff_test$SHP2_mut <- stringr::str_remove(treatment_diff_test$term,
                                                    "ClusterDEG")

treatment_diff_test$q_value <- p.adjust(p = treatment_diff_test$p_value, method = "BH")  

saveRDS(treatment_diff_test, "exclude_KOcontrols/Cluster5_vs_Q510K-T507K_full_treatment_diff_test.RDS")
treatment_diff_test <- readRDS("exclude_KOcontrols/Cluster5_vs_Q510K-T507K_full_treatment_diff_test.RDS")

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "ClusterDEG", term)) %>%
  select(Cluster = SHP2_mut, id, gene_short_name, normalized_effect, q_value)

readr::write_csv(file = "exclude_KOcontrols/Cluster5_vs_Q510K-T507K_MrVI_DEGs.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)

ggplot(treatment_diff_test %>% 
         filter(grepl(pattern = "ClusterDEG", term)) %>%
         mutate(plot_q_value = -log10(q_value)) %>%
         mutate(plot_gene = case_when(
           gene_short_name %in% "EGFR" ~ gene_short_name,
           TRUE ~ NA
         )) %>%
         mutate(plot_color = case_when(
           !is.na(plot_gene) ~ "red",
           TRUE ~ "black"
         )) %>%
         arrange(plot_color) %>%
         mutate(plot_color = factor(plot_color)),
       aes(x = normalized_effect, y = plot_q_value)) +
  # facet_wrap(~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             # alpha = 0.4,
             # stroke = 0.01,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3,
                           min.segment.length = 0.05, 
                           # force_pull = 0.4,
                           # segment.size = 0.1, box.padding = 0.1
  ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
# ggsave("SHP2_volcano_plot_KRAS_UP_highlight.png",
#        dpi = 600, height = 2.5, width = 3.5)


# Running cluster 5 against wild type not in cluster 5
colData(cds_subset)$ClusterDEG <- case_when(
  colData(cds_subset)$Cluster == 5 ~ "Cluster5",
  colData(cds_subset)$SHP2_mut %in% c("WT") ~ "Other",
  TRUE ~ NA
)

colData(cds_subset)$ClusterDEG <- factor(colData(cds_subset)$ClusterDEG, 
                                         levels = c("Other", "Cluster5"))

treatment_diff_test <- fit_models(cds = cds_subset[expressed_genes,], 
                                  model_formula_str = "~ ClusterDEG",
                                  expression_family = "quasipoisson",
                                  cores = 1, 
                                  verbose = TRUE)
treatment_diff_test <- coefficient_table(treatment_diff_test) %>% dplyr::select(-model,-model_summary)
treatment_diff_test$SHP2_mut <- stringr::str_remove(treatment_diff_test$term,
                                                    "ClusterDEG")

treatment_diff_test$q_value <- p.adjust(p = treatment_diff_test$p_value, method = "BH")  

saveRDS(treatment_diff_test, "exclude_KOcontrols/Cluster5_vsWT_diff_test.RDS")

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "ClusterDEG", term)) %>%
  select(Cluster = SHP2_mut, id, gene_short_name, normalized_effect, q_value)

readr::write_csv(file = "exclude_KOcontrols/Cluster5_vs_WT_MrVI_DEGs.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)

# ================================================================================
# Looping through clusters for DEG test (cluster n vs. all other cells) 
# ================================================================================

cds_subset <- readRDS("SHP2_cds_noKO_clustered_z-space.RDS")

treatment_diff_test.list <- list()
for (clust in 1:5) {
  colData(cds_subset)$ClusterDEG <- case_when(
    colData(cds_subset)$Cluster == clust ~ paste0("Cluster",clust),
    TRUE ~ "Other"
  )
  colData(cds_subset)$ClusterDEG <- factor(colData(cds_subset)$ClusterDEG, 
                                           levels = c("Other", paste0("Cluster",clust)))
  
  treatment_diff_test.list[[clust]] <- fit_models(cds = cds_subset[expressed_genes,], 
                                    model_formula_str = "~ ClusterDEG",
                                    expression_family = "quasipoisson",
                                    cores = 1, 
                                    verbose = TRUE)
  treatment_diff_test.list[[clust]] <- coefficient_table(treatment_diff_test.list[[clust]]) %>% 
    dplyr::select(-model,-model_summary)
  treatment_diff_test.list[[clust]]$Cluster <- stringr::str_remove(treatment_diff_test.list[[clust]]$term,
                                                                    "ClusterDEG")
  
  message("Finished processing Cluster ", clust)
  
}

treatment_diff_test <- do.call("rbind", treatment_diff_test.list)
treatment_diff_test$q_value <- p.adjust(p = treatment_diff_test$p_value, method = "BH")  

saveRDS(treatment_diff_test, "exclude_KOcontrols/Cluster_vs_Other_full_treatment_diff_test.RDS")
treatment_diff_test <- readRDS("exclude_KOcontrols/Cluster_vs_Other_full_treatment_diff_test.RDS")

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "ClusterDEG", term)) %>%
  select(Cluster, id, gene_short_name, normalized_effect, q_value)

readr::write_csv(file = "exclude_KOcontrols/Cluster_vs_Other_MrVI_DEGs.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)


library(devtools)
library(parallel)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("Q510_SHP2_EGFstim_cds_QCfiltered.RDS")

cds <- cds[,colData(cds)$SHP2_mut != "HEK293"]

colData(cds)$SHP2_mut <- factor(colData(cds)$SHP2_mut, levels = c("WT", "KO",
                                                                  "T507K", "Q510E",
                                                                  "Q510H", "Q510K",
                                                                  "Q510L", "Q510P",
                                                                  "Q510R", "Q506P", "R4A-R5A"))

colData(cds)$timepoint <- factor(colData(cds)$timepoint, levels = c("24hr", "96hr"))

colData(cds)$EGF_log10_dose <- log10(as.double(colData(cds)$EGF_dose) + 0.1)

colData(cds)$replicate <- case_when(
  colData(cds)$hash_plate == "01" ~ 1,
  colData(cds)$hash_plate == "02" ~ 2
)

colData(cds)$EGF_group <- case_when(
  colData(cds)$EGF_dose %in% c("0") ~ "No EGF",
  colData(cds)$EGF_dose %in% c("12.5", "25", "50") ~ "Low",
  colData(cds)$EGF_dose %in% c("100", "250", "500", "1000") ~ "High",
)

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= ncol(cds)*0.05,]$id

cds <- cds[,colData(cds)$SHP2_mut %in% c("WT", "Q510E", "Q510K", "Q510R")]

treatment_diff_test.list <- list()
for(time in c("24hr", "96hr")) {
  
  treatment_diff_test.list[[time]] <- fit_models(cds = cds[expressed_genes,
                                              colData(cds)$timepoint == time], 
                                    model_formula_str = "~ SHP2_mut + EGF_log10_dose + replicate",
                                    expression_family = "quasipoisson",
                                    cores = 1, 
                                    verbose = TRUE)
  treatment_diff_test.list[[time]] <- coefficient_table(treatment_diff_test.list[[time]]) %>% 
    dplyr::select(-model,-model_summary)
  treatment_diff_test.list[[time]]$SHP2_mut <- stringr::str_remove(treatment_diff_test.list[[time]]$term,
                                                      "SHP2_mut")
  treatment_diff_test.list[[time]]$timepoint <- time
  
}

treatment_diff_test <- do.call("rbind", 
                               treatment_diff_test.list)

saveRDS(treatment_diff_test, "SHP2_time_treatment_diff_test.RDS")

treatment_diff_test <- readRDS("SHP2_time_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(q_value <= 0.05)

treatment_diff_test_SHP2_plot <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  mutate(plot_color = case_when(
    abs(normalized_effect) >= 0.05 & q_value <= 0.01 ~ "red",
    TRUE ~ "black"
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T507K", "Q510E",
                                                "Q510H", "Q510K",
                                                "Q510L", "Q510P",
                                                "Q510R", "Q506P", "R4A-R5A"))) %>%
  mutate(plot_q_value = case_when(
    -log10(q_value) >= 30 ~ 30,
    TRUE ~ -log10(q_value)
  ))

ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name == "PTPN11" ~ "PTPN11",
           TRUE ~ NA
         )),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_grid(timepoint~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.1,
             show.legend = F) +
  # ggrepel::geom_label_repel(aes(label = plot_gene), 
  #                          size = 1.3,
  #                          min.segment.length = 0.05, 
  #                          max.overlaps = 2,
  #                          force_pull = 0.7,
  #                          segment.size = 0.1, box.padding = 0.1
  #                          ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  geom_hline(yintercept = -log10(0.01), linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_volcano_plot_highlight.png",
       dpi = 600, height = 2.5, width = 3)

ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name == "PTPN11" ~ "PTPN11",
           TRUE ~ NA
         )) %>%
         mutate(plot_gene_color = case_when(
           is.na(plot_gene) ~ "black",
           TRUE ~ "red"
         )) %>%
         arrange(plot_gene_color),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_grid(timepoint~SHP2_mut) +
  geom_point(aes(color = plot_gene_color),
             size = 0.1,
             show.legend = F) +
  ggrepel::geom_label_repel(aes(label = plot_gene),
                           size = 1,
                           min.segment.length = 0.05,
                           max.overlaps = 2,
                           force_pull = 0.7,
                           segment.size = 0.1, box.padding = 0.1
                           ) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  geom_hline(yintercept = -log10(0.01), linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black","red"), na.value = "black") +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_volcano_plot_PTPN11_highlight.png",
       dpi = 600, height = 2.5, width = 3)


ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           grepl("PTPN", gene_short_name) ~ gene_short_name,
           TRUE ~ NA
         )) %>%
         mutate(plot_color = case_when(
           !is.na(plot_gene) ~ "red",
           TRUE ~ "black"
         )) %>%
         arrange(plot_color) %>%
         mutate(plot_color = factor(plot_color)),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_wrap(~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3,
                           min.segment.length = 0.1, 
                           force_pull = 0.7,
                           segment.size = 0.1, box.padding = 0.1
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
ggsave("SHP2_volcano_plot_PTPN_highlight.png",
       dpi = 600, height = 2.5, width = 3.5)


treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.01) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  group_by(timepoint, SHP2_mut, direction) %>%
  dplyr::summarise(count = n()) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T507K", "Q510E",
                                                "Q510H", "Q510K",
                                                "Q510L", "Q510P",
                                                "Q510R", "Q506P", "R4A-R5A")))

ggplot(treatment_diff_test_SHP2, aes(x = SHP2_mut, y = count)) +
  facet_wrap(~timepoint, ncol = 1, scales = "free_y") +
  geom_bar(stat = "identity", 
           aes(fill = direction),
           position = "dodge") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, "Set1")[c(2,1)]) +
  xlab("SHP2 Mutant") +
  ylab("Count") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0,0,l = -0.5,0)) +
  guides(fill = guide_legend(title = "FDR < 1%\nabs(beta coeff) > 0.05")) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_mutant_significant_genes.png",
       dpi = 600, height = 4, width = 3.5)

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  # filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
  select(SHP2_mut, timepoint, id, gene_short_name, normalized_effect, q_value, p_value)

readr::write_csv(file = "SHP2_Q510_subset_vs_WT_unfiltered_mutant_diff_test_results.csv", 
          x = diff_test_for_printing %>% as.data.frame(), col_names = T)

# ================================================================================
# Correlating DEGs of mutants at 24hrs
# ================================================================================
DEG_union <- treatment_diff_test_SHP2 <- treatment_diff_test %>%
  # filter(timepoint == "24hr") %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.01) %>%
  distinct(id) %>%
  pull()

corr_df <- data.frame()
for (i in c("Q510E","Q510K","Q510R")) {
  for (j in c("Q510E","Q510K","Q510R")) {
    
    if (i == j) {
      corr_row_temp <- data.frame("SHP2_mut1" = i, "SHP2_mut2" = j, "pearson" = 1)
      corr_df <- bind_rows(corr_df, corr_row_temp)
    } else {
      treatment_diff_test_signif <- treatment_diff_test %>%
        filter(timepoint == "24hr") %>%
        filter(id %in% DEG_union) %>%
        filter(SHP2_mut %in% c(i,j)) %>%
        select(gene_id, SHP2_mut, normalized_effect) %>%
        pivot_wider(id_cols = gene_id, 
                    names_from = SHP2_mut, 
                    values_from = normalized_effect) %>%
        tibble::column_to_rownames("gene_id")
      
      corr_temp <- cor(treatment_diff_test_signif)

      corr_plot <- ggplot(treatment_diff_test_signif[c(i,j)] %>%
                            select(A = 1, B = 2) %>%
                            tibble::rownames_to_column("gene") %>%
                            mutate(diff = log10(abs(A/B))) %>%
                            arrange(desc(abs(diff))) %>%
                            dplyr::mutate(rank = row_number()) %>%
                            mutate(gene_label = case_when(
                              rank %in% c(1,2,3,4) ~ gene
                            )),
                          aes(x = A, y = B)) +
        geom_point(size = 0.1) +
        geom_smooth(method = "lm", se = FALSE, linewidth = 0.2, linetype = 2) +
        # ggrepel::geom_text_repel(aes(label = gene_label), size = 2, 
        #                          min.segment.length = 0.05, 
        #                          segment.size = 0.2, force_pull = 3) +
        annotate(geom = "text",
                 x = -0.6, y = 2.5, size = 2,
                 label = paste0("r = ", round(corr_temp[1,2], digits = 3))) +
        xlab(paste0(i)) +
        ylab(paste0(j)) +
        theme(text = element_text(size = 6),
              axis.ticks.length = unit(1,"mm"),
              axis.ticks = element_line(linewidth = 0.3)) +
        monocle3:::monocle_theme_opts()
      ggsave(plot = corr_plot,
             filename = paste0("corr_plots/",i,"_",j,"_DEG_correlation.png"),
             dpi = 600, height = 2, width = 2)
      
      corr_row_temp <- data.frame("SHP2_mut1" = i, "SHP2_mut2" = j, "pearson" = corr_temp[1,2])
      corr_df <- bind_rows(corr_df, corr_row_temp)
    }
  }
}
temp_mat <- corr_df %>%
  pivot_wider(id_cols = SHP2_mut2, 
              names_from = SHP2_mut1, 
              values_from = pearson) %>%
  tibble::column_to_rownames("SHP2_mut2") 

# hmcols <- circlize::colorRamp2(breaks = c(0,1), colors = c("white", "red"))
hm <- ComplexHeatmap::Heatmap(matrix = temp_mat %>% as.matrix(),
                              clustering_method_columns = "ward.D2",
                              clustering_method_rows = "ward.D2",
                              # col = colorRampPalette(c("blue","white","red"))(50),
                              col = viridis::plasma(n = 10),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid.text(sprintf("%.2f", temp_mat[i, j]), x, y, gp = gpar(fontsize = 5))},
                              name = "Pearson (r)",
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6),
                                                          labels_gp = gpar(fontsize = 6),
                                                          border = T,
                                                          grid_height = unit(1, "cm"),
                                                          grid_width = unit(0.4, "cm")),
                              row_title = NULL,
                              column_title_gp = gpar(fontsize = 6), 
                              # column_title_side = "bottom", column_title_rot = 45,
                              # row_title_gp = gpar(fontsize = 6),
                              # row_title_rot = 0,
                              show_column_names = T, 
                              column_names_gp = gpar(fontsize = 6),
                              column_names_rot = 45,
                              column_names_side = "bottom",
                              show_row_names = T,
                              row_names_gp = gpar(fontsize = 6),
                              cluster_rows = T,
                              cluster_columns = T,
                              rect_gp = gpar(lwd = 0.5),
                              border = T,
                              border_gp = gpar(lwd = 0.5),
                              row_dend_width = unit(0.4, "cm"),
                              column_dend_height = unit(0.4, "cm"),
                              use_raster = FALSE
)

draw(hm, merge_legends = TRUE)

png(paste0("SHP2_24hr_beta_correlation_coefficient_coeffanno.png"),
    width = 2,height = 1.5,units="in",res=1200)
ht <- draw(hm, merge_legends = TRUE)
dev.off()


temp_mat_for_printing <- corr_df %>%
  pivot_wider(id_cols = SHP2_mut2, 
              names_from = SHP2_mut1, 
              values_from = pearson) %>%
  dplyr::rename(SHP2_mut = SHP2_mut2)

write_csv(temp_mat_for_printing,
          "SHP2_Q510_select_mut_24hr_beta_coefficient_correlation_coeff_matrix.csv",
          col_names = T)


# ================================================================================
# Intersections of DEGs
# ================================================================================
treatment_diff_test_SHP2 <- treatment_diff_test %>%
  # filter(timepoint == "24hr") %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.01) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  ))

upset.list <- list()
for (time in c("24hr", "96hr")) {
  for (mut in unique(treatment_diff_test_SHP2$SHP2_mut) %>% sort()) {
    upset.list[[paste0(mut,"_",time)]] <- treatment_diff_test_SHP2 %>% 
      filter(SHP2_mut == mut & direction == "Upregulated" & timepoint == time) %>% 
      pull(gene_short_name)
  }
}


intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) > 3], pt_size = unit(1, "mm"), lwd = 0.6,
                                        comb_order = order(comb_size(intersection[comb_size(intersection) > 3])), 
                                        top_annotation = upset_top_annotation(intersection[comb_size(intersection) > 3],
                                                                              # ylim = c(0, 500),
                                                                              add_numbers = T,
                                                                              numbers_rot = 0,
                                                                              numbers_gp = gpar(fontsize = 6),
                                                                              height = unit(1, "cm"),
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              annotation_name_gp = gpar(fontsize = 6),
                                                                              show_annotation_name = T),
                                        right_annotation = upset_right_annotation(intersection[comb_size(intersection) > 3],
                                                                                  add_numbers = FALSE,
                                                                                  width = unit(1, "cm"),
                                                                                  axis_param = list(gp = gpar(fontsize = 6)),
                                                                                  annotation_name_gp = gpar(fontsize = 6),
                                                                                  show_annotation_name = T),
                                        row_names_gp = gpar(fontsize = 6, fontface = "bold"), 
                                        column_title = "SHP2mut DEG Intersection (upregulated, FDR < 0.01)",
                                        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

draw(intersect_plot)

png("SHP2_mut_DEG_upreg_intersection.png",width=4.5,height=2,units="in",res=1200)
draw(intersect_plot)
dev.off()

intersection <- intersection[comb_size(intersection) >= 3]
gene_sets <- sapply(comb_name(intersection), function(nm) extract_comb(intersection, nm))
test <- set_name(intersection)
test_names <- data.frame(x = names(gene_sets)) %>%
  separate(x, into = as.character(c(1:6)), sep = c(1,2,3,4,5,6)) %>%
  mutate(`1` = ifelse(`1` == 1, test[1],""),
         `2` = ifelse(`2` == 1, test[2],""),
         `3` = ifelse(`3` == 1, test[3],""),
         `4` = ifelse(`4` == 1, test[4],""),
         `5` = ifelse(`5` == 1, test[5],""),
         `6` = ifelse(`6` == 1, test[6],"")) %>%
  unite(col = "name", sep = " ", 1:6)

clean_mutations <- function(x) {
  x <- gsub("^\\s+|\\s+$", "", x)     # Remove leading and trailing spaces
  x <- gsub("\\s+", " ", x)           # Replace multiple spaces with a single space
  x <- gsub(" ", "_", x)              # Replace spaces with underscores
  return(x)
}

# Apply the function to the name column
name <- sapply(test_names, clean_mutations)

# Convert back to data frame if necessary
names(gene_sets) <- name

write_tsv(x = as.data.frame(gene_sets["Q510K_24hr_Q510R_24hr"]), 
          file = "SHP2_Q510_K_and_R_DEG_intersection_vs_WT.txt")

# Grab DEG test information for each element in list
gene_sets_final <- list()
for (new_element in names(gene_sets)) {
  new_element_split <- stringr::str_split_1(string = new_element, pattern = "_")
  gene_set_df_temp <- treatment_diff_test_SHP2 %>%
    filter(gene_short_name %in% gene_sets[[new_element]]) %>%
    filter(SHP2_mut %in% new_element_split) %>%
    select(-id, -status, -estimate, -std_err, -test_val, -p_value, -model_component) %>%
    as.data.frame()
  gene_sets_final[[new_element]] <- gene_set_df_temp
}

# ================================================================================
# Heatmaps of gene intersections
# ================================================================================

gene_set_of_interest <- gene_sets["Q510K_24hr_Q510R_24hr"] %>% unlist()
# gene_set_of_interest <- gene_sets["Q510E_24hr"] %>% unlist()

cds_subset <- cds[rowData(cds)$gene_short_name %in% gene_set_of_interest,
                  colData(cds)$timepoint == "24hr" &
                    # colData(cds)$EGF_group == "No EGF"]
                    colData(cds)$SHP2_mut != "WT"]

cell_group_df <- colData(cds_subset) %>% as.data.frame() %>%
  select(cell = cell_ID,
         SHP2_mut,
         EGF_group) %>%
  unite("Group", SHP2_mut, EGF_group, sep = "_")
agg_mat <- aggregate_gene_expression(cds_subset, gene_group_df = NULL, 
                                     cell_group_df, cell_agg_fun = "mean")
row_names_temp <- data.frame("temp" = row.names(agg_mat)) %>%
  left_join(rowData(cds) %>% as.data.frame() %>% select(id, gene_short_name),
            by = c("temp" = "id")) %>%
  select(gene_short_name) %>%
  pull()
row.names(agg_mat) <- row_names_temp
colnames(agg_mat) <- stringr::str_c(colnames(agg_mat))

col_anno_df <- data.frame("temp" = colnames(agg_mat)) %>%
  separate(col = temp, into = c("SHP2", "EGF Dose"), sep = "_", remove = F) %>%
  tibble::column_to_rownames("temp") %>%
  mutate(`EGF Dose` = factor(`EGF Dose`, levels = c("No EGF", "Low", "High")))

col_anno <- rowAnnotation(df = col_anno_df,
                              col = list(`SHP2` = c("Q510E" = "#F8766D",
                                                    "Q510K" = "#00BFC4",
                                                    "Q510R" = "purple"),
                                         `EGF Dose` = c("No EGF" = "grey90",
                                                        "Low" = viridis::magma(n = 4)[2],
                                                        "High" = viridis::magma(n = 4)[4])), 
                              show_legend = T,
                              annotation_legend_param = list(`SHP2` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                                           labels_gp = gpar(fontsize = 6)),
                                                             `EGF Dose` = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                                               labels_gp = gpar(fontsize = 6))), 
                              annotation_name_gp = gpar(fontsize = 6)
)

hmcol <- circlize::colorRamp2(c(2,1,0,-1,-2), RColorBrewer::brewer.pal(n = 5, "RdYlBu"))
hm <- ComplexHeatmap::Heatmap(matrix = agg_mat %>% t(),
                              name = "Z-Scored\nMean Epxression",
                              col = hmcol,
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 6, fontface = "bold"),
                                                          labels_gp = gpar(fontsize = 6),
                                                          border = TRUE),
                              clustering_method_columns = "ward.D2", 
                              clustering_method_rows = "ward.D2",
                              # right_annotation = col_anno,
                              row_title = "SHP2_EGF Group",
                              column_title_side = "bottom",
                              column_title_gp = gpar(fontsize = 8, fontface = "bold"),
                              column_title = "Gene",
                              row_title_side = "left",
                              row_title_gp = gpar(fontsize = 8, fontface = "bold"),
                              show_column_names = T,
                              column_names_rot = 45,
                              column_names_gp = gpar(fontsize = 6),
                              show_row_names = T,
                              row_names_gp = gpar(fontsize = 6),
                              row_names_side = "left",
                              cluster_rows = T,
                              cluster_columns = T,
                              row_split = 3,
                              column_split = 4,
                              row_dend_width = unit(0.2,"cm"),
                              column_dend_height = unit(0.4,"cm"),
                              use_raster = FALSE
)

draw(hm, legend_grouping = "original")

pdf(file = "SHP2_mut_DEG_upreg_Q510K_Q510R_intersection_expression.pdf",
    width = 6.5, height = 2.5, compress = F, bg = "white")
draw(hm, legend_grouping = "original")
dev.off()


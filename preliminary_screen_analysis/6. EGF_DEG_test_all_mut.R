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

cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

colData(cds)$SHP2_mut <- factor(colData(cds)$SHP2_mut, levels = c("KO", "WT",
                                                                  "T42A", "T52S",
                                                                  "E76K", "R138Q",
                                                                  "E139D", "Y279C",
                                                                  "T468M", "T507K",
                                                                  "Q510E", "Q510K"))

colData(cds)$timepoint <- factor(colData(cds)$timepoint, levels = c("24hr", "96hr"))

colData(cds)$EGF_log10_dose <- log10(as.double(colData(cds)$EGF_dose) + 0.1)

colData(cds)$replicate <- case_when(
  colData(cds)$hash_plate == "01" ~ 1,
  colData(cds)$hash_plate == "02" ~ 2
)

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.01,]$id

treatment_diff_test.list <- list()
for (mut in c("KO", "WT",
              "T42A", "T52S",
              "E76K", "R138Q",
              "E139D", "Y279C",
              "T468M", "T507K",
              "Q510E", "Q510K")) {
  treatment_diff_test.list[[mut]] <- fit_models(cds = cds[expressed_genes,
                                                   colData(cds)$SHP2_mut == mut], 
                                    model_formula_str = "~ timepoint + EGF_log10_dose + replicate",
                                    expression_family = "quasipoisson",
                                    cores = 1, 
                                    verbose = TRUE)
  treatment_diff_test.list[[mut]] <- coefficient_table(treatment_diff_test.list[[mut]]) %>% 
    dplyr::select(-model,-model_summary)
  treatment_diff_test.list[[mut]] <- treatment_diff_test.list[[mut]] %>%
    mutate(SHP2_mut = mut)
  message(paste0("Finished processing SHP2_", mut))
}

treatment_diff_test <- do.call("rbind", 
                               treatment_diff_test.list)

# saveRDS(treatment_diff_test, "SHP2_mut_EGF_treatment_diff_test.RDS")
# treatment_diff_test <- readRDS("SHP2_mut_EGF_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(q_value <= 0.05)

treatment_diff_test_SHP2_plot <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  mutate(plot_color = case_when(
    abs(normalized_effect) >= 0.05 & q_value <= 0.05 ~ "red",
    TRUE ~ "black"
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K",
                                                "Q510E", "Q510K"))) %>%
  mutate(plot_q_value = case_when(
    -log10(q_value) >= 10 ~ 10,
    TRUE ~ -log10(q_value)
  )) %>%
  mutate(plot_normalized_effect = case_when(
    normalized_effect >= 0.6 ~ 0.6,
    normalized_effect <= -0.6 ~ -0.6,
    TRUE ~ normalized_effect
  ))


ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name == "PTPN11" ~ "PTPN11",
           TRUE ~ NA
         )) %>%
         mutate(SHP2_mut = case_when(
           SHP2_mut == "Q510E" ~ "Transfection\nControl",
           TRUE ~ SHP2_mut
         )) %>%
         mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "Transfection\nControl",
                                                       "WT", "T42A", "T52S",
                                                       "E76K", "R138Q",
                                                       "E139D", "Y279C",
                                                       "T468M", "T507K","Q510K"))),
       aes(x = plot_normalized_effect, y = plot_q_value)) +
  facet_wrap(~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  geom_vline(xintercept = 0.05, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.05, linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black", "red")) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("EGF_volcano_plot_highlight.png",
       dpi = 600, height = 3, width = 3)

ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name == "EGFR" ~ gene_short_name,
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
                           force_pull = 0.5,
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
ggsave("SHP2_volcano_plot_EGFR_highlight.png",
       dpi = 600, height = 3, width = 3)

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  group_by(SHP2_mut, direction) %>%
  dplyr::summarise(count = n()) %>%
  mutate(SHP2_mut = case_when(
    SHP2_mut == "Q510E" ~ "Transfection\nControl",
    TRUE ~ SHP2_mut
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "Transfection\nControl",
                                                "WT", "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K","Q510K")))

ggplot(treatment_diff_test_SHP2, aes(x = SHP2_mut, y = count)) +
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
  guides(fill = guide_legend(title = "FDR < 5%\nabs(beta coeff) > 0.05")) +
  monocle3:::monocle_theme_opts()
ggsave("EGF_mutant_significant_genes.png",
       dpi = 600, height = 1.5, width = 3.5)

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  # filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
  select(SHP2_mut, id, gene_short_name, normalized_effect, q_value, p_value)

readr::write_csv(file = "SHP2_EGF_unfiltered_EGFdose_diff_test_results.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "timepoint96", term)) %>%
  # filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
  select(SHP2_mut, term, id, gene_short_name, normalized_effect, q_value, p_value)

readr::write_csv(file = "SHP2_EGF_unfiltered_timepoint_diff_test_results.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)


# ==============================================================================
# Analyzing intersections of DEGs between SHP2 mutants
# ==============================================================================

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  mutate(SHP2_mut = case_when(
    SHP2_mut == "Q510E" ~ "Transfection\nControl",
    TRUE ~ SHP2_mut
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "Transfection\nControl",
                                                "WT", "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K","Q510K")))

DEG_union <- treatment_diff_test_SHP2 %>%
  distinct(id) %>%
  pull()

# correlation matrix
corr_df <- data.frame()
  for (i in c("KO", "Q510E",
              "WT", "T42A", "T52S",
              "E76K", "R138Q",
              "E139D", "Y279C",
              "T468M", "T507K","Q510K")) {
    for (j in c("KO", "Q510E",
                "WT", "T42A", "T52S",
                "E76K", "R138Q",
                "E139D", "Y279C",
                "T468M", "T507K","Q510K")) {
      
      if (i == j) {
        corr_row_temp <- data.frame("SHP2_mut1" = i, "SHP2_mut2" = j, "pearson" = 1)
        corr_df <- bind_rows(corr_df, corr_row_temp)
      } else {
        treatment_diff_test_signif <- treatment_diff_test %>%
          filter(term %in% c("EGF_log10_dose")) %>%
          # filter(id %in% DEG_union) %>%
          filter(SHP2_mut %in% c(i,j)) %>%
          select(gene_id, SHP2_mut, normalized_effect) %>%
          pivot_wider(id_cols = gene_id, names_from = SHP2_mut, values_from = normalized_effect) %>%
          tibble::column_to_rownames("gene_id")
        
        corr_temp <- cor(treatment_diff_test_signif)
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
                                col = colorRampPalette(c("white","red"))(50),
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
                                use_raster = FALSE
  )
  
  draw(hm, merge_legends = TRUE)
  
  png(paste0("beta_correlation_coefficient_",
             PDCL,
             ".png"),
      width = 3,height = 2.5,units="in",res=1200)
  ht <- draw(hm, merge_legends = TRUE)
  dev.off()



# ==============================================================================
# GSEA of SHP2mut EGF DEGs
# ==============================================================================
treatment_diff_test <- readRDS("SHP2_mut_EGF_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")


library(fgsea)
# https://stephenturner.github.io/deseq-to-fgsea/
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

pathways.hallmark <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/h.all.v6.0.symbols.gmt")
pathways.reactome <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
# pathways.gobp <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
pathways.oncogenic <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c6.all.v2023.2.Hs.symbols.gmt")

# ranks must look like a named vector of values
##         A1BG     A1BG-AS1         A1CF          A2M      A2M-AS1 
##  0.679946437 -1.793291412 -0.126192849 -1.259539478  0.875346116 
# data(exampleRanks)

# DEG_union <- treatment_diff_test_SHP2 %>%
#   distinct(id) %>%
#   pull()

library(gridExtra)
fgseaRes_list <- list()
for (group in c("KO", "WT",
                "T42A", "T52S",
                "E76K", "R138Q",
                "E139D", "Y279C",
                "T468M", "T507K",
                "Q510E", "Q510K")) {
  ranks <- treatment_diff_test %>%
    filter(grepl(pattern = "EGF", term)) %>%
    filter(SHP2_mut == group) %>%
    # filter(id %in% DEG_union) %>%
    arrange(desc(normalized_effect)) %>%
    select(gene_short_name, normalized_effect)
  
  rank_vector <- ranks$normalized_effect
  names(rank_vector) <- ranks$gene_short_name
  
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=rank_vector)
  
  fgseaRes_list[[group]] <- fgseaRes %>% mutate(SHP2_mut = group)
  
}

fgseaRes_final <- do.call("rbind",fgseaRes_list)
fgseaRes_final$q_value <- p.adjust(fgseaRes_final$pval, method = "BH")

fgseaRes_plot <- fgseaRes_final %>%
  filter(q_value <= 0.01)

ggplot(fgseaRes_plot, aes(x = NES, y = pathway)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = SHP2_mut)) +
  # geom_point(aes(color = SHP2_mut)) +
  monocle3:::monocle_theme_opts()



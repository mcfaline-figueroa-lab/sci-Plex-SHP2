library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
library(gghalves)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)
source("calculate_aggreg_expression.R")

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

p <- plot_percent_cells_positive(cds[rowData(cds)$gene_short_name %in% c("PTPN11"),
                                colData(cds)$SHP2_mut %in% c("KO", "WT") & 
                                  colData(cds)$timepoint == "24hr"],
                            group_cells_by = "SHP2_mut", min_expr = 0) +
  ylab("% Positive Cells") +
  theme(text = element_text(size = 9, color = "black"),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(1,"mm"),
        axis.title.x = element_blank()) +
  # scale_y_continuous(breaks = seq(0,25,5)) +
  guides(fill = "none")
p
ggsave(plot = p, "SHP2_percent_cells_positive.png",dpi = 600, height = 1.5, width = 1.5)


colData(cds)$PTPN11 <- calculate_aggregate_expression_score(cds,signature_genes = "PTPN11")
plot_data <- colData(cds) %>% 
  as.data.frame() %>% 
  filter(timepoint == "24hr") %>%
  filter(SHP2_mut %in% c("KO", "WT"))

plot_data_mean <- plot_data %>%
  dplyr::group_by(SHP2_mut) %>%
  dplyr::summarise(mean_oxphos = mean(log2(PTPN11 + 1)), median_oxphos = median(PTPN11)) %>%
  mutate(x_axis = case_when(
    SHP2_mut == "WT" ~ 0.05,
    TRUE ~ -0.05
  ))

ggplot(plot_data,
       aes(y = log2(PTPN11 + 1))) +
  geom_half_violin(aes(split = SHP2_mut, fill = SHP2_mut),
                   show.legend = F, 
                   nudge = -.19,
                   scale = "width",
                   linewidth = 0.2) +
  geom_point(data = plot_data_mean, aes(x = x_axis, y = mean_oxphos), size = 0.2) +
  scale_x_continuous(breaks = c(-0.14, 0.14), labels = c("KO", "WT")) +
  # stat_summary(fun.y=mean, geom="point", color="black", stroke = 0, size = 0.8) +
  # facet_grid(timepoint~EGF_group, scales = "free") +
  # facet_wrap(~EGF_group, ncol = 3) +
  monocle3:::monocle_theme_opts() +
  ylab("Oxidative Phosphorylation\nAggregate Score") +
  theme(text = element_text(size = 8),
        # axis.ticks.length.y = unit(0.75,"mm"),
        # axis.ticks.y = element_line(linewidth = 0.2),
        # axis.ticks.x = element_blank(),
        # axis.text.x = element_blank(),
        axis.ticks.length = unit(0.75,"mm"),
        axis.ticks = element_line(linewidth = 0.2),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black")
  )

treatment_diff_test <- readRDS("../SHP2_mut_EGF_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(q_value <= 0.05)

treatment_diff_test_SHP2_plot <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  mutate(plot_color = case_when(
    (normalized_effect) >= 0.1 & q_value <= 0.05 ~ "red",
    (normalized_effect) <= -0.1 & q_value <= 0.05 ~ "blue",
    TRUE ~ "black"
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K",
                                                "Q510E", "Q510K"))) %>%
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
  facet_wrap(~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.1,
             show.legend = F) +
  geom_vline(xintercept = 0.1, linetype = 2, linewidth = 0.1) +
  geom_vline(xintercept = -0.1, linetype = 2, linewidth = 0.1) +
  scale_color_manual(values = c("black",RColorBrewer::brewer.pal(n = 3, "Set1")[c(2,1)])) +
  xlab("Normalized Effect") +
  ylab("-Log10 FDR") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1)) +
  monocle3:::monocle_theme_opts()
ggsave("EGF_volcano_plot_highlight.png",
       dpi = 600, height = 1.5, width = 2)

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
ggsave("SHP2_volcano_plot_EGFR_highlight.png",
       dpi = 600, height = 4, width = 4)

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  # filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  dplyr::group_by(SHP2_mut, direction) %>%
  dplyr::summarise(count = n()) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K",
                                                "Q510E", "Q510K")))

write_csv(x = treatment_diff_test_SHP2, file = "SHP2_EGF_filtered_EGFdose_diff_test_counts.csv")

ggplot(treatment_diff_test_SHP2, aes(x = SHP2_mut, y = count)) +
  geom_bar(stat = "identity", 
           aes(fill = direction),
           position = "dodge", width = 0.8) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, "Set1")[c(2,1)]) +
  xlab("SHP2 Mutant") +
  ylab("Count") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0,0,t = -8,0),
        legend.position = "bottom") +
  guides(fill = guide_legend(title = "FDR < 5%\nabs(beta coeff) > 0.1", direction = "horizontal", ncol = 1)) +
  monocle3:::monocle_theme_opts()
ggsave("EGF_mutant_significant_genes.png",
       dpi = 600, height = 1.5, width = 1.75)

diff_test_for_printing <- treatment_diff_test %>%
  filter(grepl(pattern = "EGF", term)) %>%
  # filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
  select(SHP2_mut, id, gene_short_name, normalized_effect, q_value, p_value)

readr::write_csv(file = "SHP2_unfiltered_mutant_EGF_diff_test.csv", 
                 x = diff_test_for_printing %>% as.data.frame(), col_names = T)

cutoff_df <- data.frame()
for (cutoff in c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)) {
  treatment_diff_test_SHP2 <- treatment_diff_test %>%
    filter(SHP2_mut %in% c("KO", "WT")) %>%
    filter(grepl(pattern = "EGF", term)) %>%
    filter(abs(normalized_effect) >= cutoff & q_value <= 0.05) %>%
    mutate(direction = case_when(
      normalized_effect >= 0 ~ "Upregulated",
      TRUE ~ "Downregulated"
    )) %>%
    group_by(SHP2_mut) %>%
    dplyr::summarise(count = n()) %>%
    mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT"))) %>%
    mutate(cutoff = cutoff)
  cutoff_df <- bind_rows(cutoff_df, treatment_diff_test_SHP2)
}

write_csv(x = cutoff_df %>%
            select(SHP2_mut, norm_effect_cutoff = cutoff, EGF_DEG_count = count), 
          file = "SHP2_WT_vs_KO_EGF_DEG_FDR0.05_norm_effect_cutoff.csv",
          col_names = T)

ggplot(cutoff_df %>% dplyr::rename(SHP2 = SHP2_mut), 
       aes(x = cutoff, y = count, group = SHP2)) +
  geom_point(aes(color = SHP2)) +
  geom_line(aes(color = SHP2)) +
  xlab("abs(Normalized Effect) Cut-off") +
  ylab("DEG Count (FDR < 5%)") +
  theme(text = element_text(size = 9),
        legend.margin = margin(l = -10)) +
  monocle3:::monocle_theme_opts()

ggsave("EGF_DEGs_WT_vs_KO_normeffect_cutoff_line.png", 
       dpi = 600, height = 2, width = 2.75)

# ==============================================================================
# Analyzing intersections of DEGs between SHP2 mutants
# ==============================================================================

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  ))

upset.list <- list()
for (mut in unique(treatment_diff_test_SHP2$SHP2_mut) %>% sort()) {
  upset.list[[paste0(mut,"_up")]] <- treatment_diff_test_SHP2 %>% 
    filter(SHP2_mut == mut & direction == "Upregulated") %>% 
    pull(id)
  upset.list[[paste0(mut,"_down")]] <- treatment_diff_test_SHP2 %>% 
    filter(SHP2_mut == mut & direction == "Downregulated") %>% 
    pull(id)
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection, pt_size = unit(2, "mm"), lwd = 0.9,
                                        comb_order = order(comb_size(intersection)), 
                                        set_order = c("WT_up",
                                                      "WT_down",
                                                      "KO_up",
                                                      "KO_down"),
                                        top_annotation = upset_top_annotation(intersection,
                                                                              # ylim = c(0, 500),
                                                                              add_numbers = T,
                                                                              numbers_rot = 0,
                                                                              numbers_gp = gpar(fontsize = 6),
                                                                              height = unit(1, "cm"),
                                                                              axis_param = list(gp = gpar(fontsize = 6)),
                                                                              annotation_name_gp = gpar(fontsize = 6),
                                                                              show_annotation_name = T),
                                        right_annotation = upset_right_annotation(intersection,
                                                                                  add_numbers = FALSE,
                                                                                  width = unit(1, "cm"),
                                                                                  axis_param = list(gp = gpar(fontsize = 6)),
                                                                                  annotation_name_gp = gpar(fontsize = 6),
                                                                                  show_annotation_name = T),
                                        row_names_gp = gpar(fontsize = 6, fontface = "bold"), 
                                        column_title = "SHP2mut DEG Intersection (FDR < 0.05)",
                                        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

png("EGF_DEG_BOTHdirections_intersection.png",width=3.5,height=1.5,units="in",res=1200)
draw(intersect_plot)
dev.off()

# correlation of DEG betas
DEG_space <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(q_value <= 0.05) %>%
  distinct(id) %>%
  mutate(signif = "Yes")

corr_plot <- treatment_diff_test %>%
  filter(SHP2_mut %in% c("KO", "WT")) %>%
  filter(grepl(pattern = "EGF", term)) %>%
  filter(id %in% DEG_space$id) %>%
  select(id, gene_short_name, normalized_effect, SHP2_mut) %>%
  pivot_wider(id_cols = c(gene_short_name, id), values_from = normalized_effect, names_from = SHP2_mut)

pearson_cor <- cor(x = corr_plot %>% select(KO, WT), method = "pearson")

ggplot(corr_plot %>% left_join(DEG_space), 
       aes(x = WT, y = KO)) +
  # geom_point(size = 0.05, aes(color = signif)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", size = 0.5, linetype = 1, se = F) +
  geom_hline(size = 0.35, yintercept = 0, linetype = 3) +
  geom_vline(size = 0.35, xintercept = 0, linetype = 3) +
  annotate(geom = "text", x = -0.3, y = 0.28, 
           label = paste0("r = ", round(pearson_cor["WT","KO"], 3)),
           size = 2) +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(1, "mm"),
        axis.ticks = element_line(linewidth = 0.5)) +
  monocle3:::monocle_theme_opts()
ggsave("EGF_dose_WT_vs_KO_deg_corr.png",
       dpi = 600, height = 1.5, width = 1.5)


colData(cds)$EGF_log10_dose <- as.factor(round(log10(as.double(colData(cds)$EGF_dose) + 1), 2))

p1 <- plot_percent_cells_positive(cds[rowData(cds)$gene_short_name == "RHOD", 
                                colData(cds)$SHP2_mut %in% c("KO")], 
                  group_cells_by = "EGF_log10_dose") +
  # scale_y_continuous(breaks = seq(0,15,3), limits = c(0,17)) +
  scale_y_continuous(breaks = seq(0,8,2), limits = c(0,9)) +
  guides(fill = "none")

p2 <- plot_percent_cells_positive(cds[rowData(cds)$gene_short_name == "RHOD", 
                                      colData(cds)$SHP2_mut %in% c("WT")], 
                                  group_cells_by = "EGF_log10_dose") +
  # scale_y_continuous(breaks = seq(0,15,3), limits = c(0,17)) +
  scale_y_continuous(breaks = seq(0,8,2), limits = c(0,9)) +
  guides(fill = "none")

cowplot::plot_grid(plotlist = list(p1,p2), nrow = 1)

# ==============================================================================
# GSEA of SHP2mut EGF DEGs
# ==============================================================================
treatment_diff_test <- readRDS("../SHP2_mut_EGF_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")


library(fgsea)
# https://stephenturner.github.io/deseq-to-fgsea/
# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

pathways.hallmark <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/h.all.v6.0.symbols.gmt")
pathways.reactome <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c2.cp.reactome.v2024.1.Hs.symbols.gmt")
pathways.gobp <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt")
pathways.oncogenic <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c6.all.v2023.2.Hs.symbols.gmt")
pathways.tft <- gmtPathways("/Users/rossgiglio/Documents/McFaline-Figueroa_Lab/GSEA_helper_functions/c3.tft.v2024.1.Hs.symbols.gmt")

# ranks must look like a named vector of values
##         A1BG     A1BG-AS1         A1CF          A2M      A2M-AS1 
##  0.679946437 -1.793291412 -0.126192849 -1.259539478  0.875346116 
# data(exampleRanks)

# DEG_union <- treatment_diff_test_SHP2 %>%
#   distinct(id) %>%
#   pull()

library(gridExtra)
fgseaRes_list <- list()
for (group in c("KO", "WT")) {
  ranks <- treatment_diff_test %>%
    filter(grepl(pattern = "EGF", term)) %>%
    filter(SHP2_mut == group) %>%
    # filter(id %in% DEG_union) %>%
    arrange(desc(normalized_effect)) %>%
    select(gene_short_name, normalized_effect) %>%
    distinct(gene_short_name, .keep_all = T)
  
  rank_vector <- ranks$normalized_effect
  names(rank_vector) <- ranks$gene_short_name
  
  fgseaRes <- fgsea(pathways=pathways.reactome, stats=rank_vector)
  fgseaRes_list[[group]] <- fgseaRes %>% mutate(SHP2_mut = group)
  
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=rank_vector)
  fgseaRes_list[[group]] <- bind_rows(fgseaRes_list[[group]], fgseaRes %>% mutate(SHP2_mut = group))
  
}

fgseaRes_final <- do.call("rbind",fgseaRes_list)
fgseaRes_final$q_value <- p.adjust(fgseaRes_final$pval, method = "BH")

fgseaRes_final <- fgseaRes_final %>%
  mutate(signif = ifelse(q_value <= 0.05, "yes", "no"))

write_csv(fgseaRes_final %>%
            select(SHP2_mut, pathway, NES, q_value, signif), 
          "SHP2_EGFdose_DEGs_WT_KO_GSEA.csv")

fgseaRes_plot <- fgseaRes_final %>%
  filter(q_value <= 0.01)

ggplot(fgseaRes_final %>% filter(pathway %in% fgseaRes_plot$pathway), 
       aes(x = NES, y = pathway, fill = SHP2_mut)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = SHP2_mut)) +
  monocle3:::monocle_theme_opts()

ggplot(fgseaRes_final %>% 
         filter(pathway %in% c("HALLMARK_MTORC1_SIGNALING",
                               "HALLMARK_MYC_TARGETS_V1",
                               "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                               "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
                               "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
                               "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION")), 
       aes(x = NES, y = pathway, fill = SHP2_mut)) +
  geom_bar(stat = "identity", position = "dodge", aes(fill = SHP2_mut)) +
  geom_text(aes(label = signif(q_value, 2), group = SHP2_mut, hjust = ifelse(NES > 0, 1.1, -0.1)), 
            position = position_dodge(width = .9),
            size = 3) +
  xlab("Normalized Effect Score (NES)") +
  theme(axis.title.y = element_blank(),
        text = element_text(size = 10),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.margin = margin(t = -7)) +
  guides(fill = guide_legend(title = "SHP2", override.aes = list(size = 2), ncol = 2)) +
  monocle3:::monocle_theme_opts()
ggsave("GSEA_selected_EGF_DEGs_WT_vs_KO.png",
       dpi = 600, height = 2.5, width = 6.5)

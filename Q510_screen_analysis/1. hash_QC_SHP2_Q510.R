library(monocle3)
library(tidyverse)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("raw_data/cds_precell_prehash.RDS")
hash <- read_tsv("raw_data/hashTable.out", col_names = c("sample", "cell_ID", "oligo", "hash_umis"))

cds_col <- colData(cds) %>% as.data.frame() %>%
  mutate(cell_ID=Cell) %>%
  tidyr::separate(col = Cell, into = c('P7', 'P5', NA, NA, 'n1', NA, NA, 'n2'), sep = "_") %>%
  mutate(RT = paste('RT_BC_', n1, sep = '')) %>%
  mutate(Lig = paste('Lig_BC_', n2, sep = '')) %>%
  select(-n1, -n2)

RT_map <- read_csv("raw_data/CellLine_RT_Map.csv", col_names = TRUE)

cds_col <- left_join(cds_col, RT_map, by = c("RT" = "RT_barcode"))

hashTable_summary <- hash %>%
  group_by(cell_ID) %>%
  dplyr::mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  dplyr::arrange(desc(proportion)) %>%
  dplyr::mutate(rank = row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  dplyr::mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, proportion[rank == 1]/proportion[rank == 2], 10)) %>%
  dplyr::filter(rank == 1)

cds_col <- left_join(cds_col, hashTable_summary, by = c("sample" = "sample", "cell_ID" = "cell_ID"))

cds_col_temp <- cds_col %>%
  separate(oligo, into = c("SHP2_mut", "EGF_dose", "hash"), sep = "_", convert = T) %>%
  separate(hash, sep = 2, into = c("hash_plate", "hash_well")) %>%
  mutate(log10_dose = log10(EGF_dose)) %>%
  dplyr::rename(n_umi = n.umi)

test <- cds_col_temp %>% group_by(SHP2_mut, EGF_dose, timepoint) %>% summarise(count = n())

colData(cds) <- DataFrame(cds_col_temp)
colData(cds)$Cell <- colData(cds)$cell_ID
row.names(colData(cds)) <- colData(cds)$cell_ID

colData(cds)$log10_umi <- colData(cds)$n_umi %>% log10()
mt_genes <- rowData(cds) %>% as.data.frame() %>%
  filter(grepl("MT-",gene_short_name) == TRUE)
mt_genes_id <- rownames(mt_genes)
colData(cds)$percent_mito <- 100 * (colSums(exprs(cds)[mt_genes_id,])/colSums(exprs(cds)))

ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  # geom_point(size = 0.5, stroke = 0, aes(color = timepoint)) +
  geom_point(size = 0.5, stroke = 0) +
  geom_hline(linewidth = 0.2, yintercept = log10(300), linetype = 2) +
  monocle3:::monocle_theme_opts() +
  # scale_color_manual(name = "Timepoint",
  #                    breaks = c("24hr", "96hr"),
  #                    values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)]) +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  annotate(geom = "text", 
           label = "cells post-UMI cutoff: 107,726",
           x = 1.5, y = 0.5,
           size = 1) +
  annotate(geom = "text", 
           label = "cutoff = 300",
           x = 1, y = log10(450),
           size = 1) +
  # guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.margin = unit(c(-20,0,0,0), "mm"),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1))
ggsave("QC/Kneeplot_300umi_allcells.png", dpi = 600, height = 1.25, width = 1.5)

write_csv(x = colData(cds) %>% 
            as.data.frame() %>% 
            arrange(desc(n_umi)) %>% 
            mutate(cell_rank = dplyr::row_number()) %>%
            mutate(log10_umi = log10(n_umi),
                   log10_cell_rank = log10(cell_rank)) %>%
            select(log10_umi, log10_cell_rank),
          "raw_data/SHP2_Q510_screen_unfiltered_for_kneeplot.csv")

p <- ggplot(colData(cds) %>% 
         as.data.frame() %>% 
           mutate(facet_option = case_when(
             SHP2_mut == "HEK293" ~ "HEK293",
             TRUE ~ "Edited"
           )) %>%
           filter(facet_option == "HEK293") %>%
         arrange(desc(n_umi)) %>% 
         mutate(cell_rank = dplyr::row_number()), 
       aes(x= log10(cell_rank), y = log10(n_umi))) +
  # geom_point(size = 0.5, stroke = 0, aes(color = timepoint)) +
  geom_point(size = 0.4, stroke = 0,
             aes(color = facet_option)) +
  geom_hline(linewidth = 0.2, yintercept = log10(300), linetype = 2) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "Background",
                     breaks = c("HEK293", "Edited"),
                     values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)]) +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  # annotate(geom = "text", 
  #          label = "cells post-UMI cutoff: 107,726",
  #          x = 1.5, y = 0.5,
  #          size = 1) +
  # annotate(geom = "text", 
  #          label = "cutoff = 300",
  #          x = 1, y = log10(450),
  #          size = 1) +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.margin = unit(c(-20,0,0,0), "mm"),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1))
ggsave(plot = p, "QC/Kneeplot_300umi_HEK293.png", dpi = 600, height = 1.5, width = 2)


g <- ggplot(colData(cds) %>% as.data.frame() %>% arrange(desc(n_umi)) %>% mutate(cell_rank = dplyr::row_number()), aes(x= log10(cell_rank), y = log10(n_umi))) +
  # geom_point(size = 0.5, stroke = 0, aes(color = timepoint)) +
  geom_point(size = 0.5, stroke = 0, aes(color = P7P5)) +
  geom_hline(linewidth = 0.2, yintercept = log10(350)) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "P7 + P5",
                     breaks = c("B02", "E06", "G07"),
                     values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2,3)],
                     labels = c("B + 02",
                                "E + 06",
                                "G + 07")) +
  scale_y_continuous(breaks = c(1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  # guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.margin = unit(c(-20,0,0,0), "mm"),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1))
ggsave(plot = g, "QC/Kneeplot_350umi_allcells.png", dpi = 600, height = 1.5, width = 2)


g <- ggplot(colData(cds) %>% as.data.frame() %>% filter(hash_umis >= 5)) +
  ggridges::geom_density_ridges(aes(x = log10(n_umi), y = P7P5, fill = P7P5), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 UMI") +
  ylab("Density") +
  scale_x_continuous(breaks = seq(0,5,1)) +
  scale_fill_manual(name = "P7 + P5",
                     breaks = c("B02", "E06", "G07"),
                     values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2,3)],
                     labels = c("B + 02",
                                "E + 06",
                                "G + 07"))
ggsave(plot = g, "QC/umi_density_300umi_hashfilt5.png", dpi = 600, height = 2, width = 2.5)

saveRDS(cds, "SHP2mut_EGFstim_Q510_cds_unfiltered.RDS")

dim(cds[,colData(cds)$n_umi >= 100 & !is.na(colData(cds)$hash_umis) & colData(cds)$hash_umis >= 3 & colData(cds)$top_to_second_best_ratio >= 2])
dim(cds[,colData(cds)$n_umi >= 300 & !is.na(colData(cds)$hash_umis) & colData(cds)$hash_umis >= 5 & colData(cds)$top_to_second_best_ratio >= 2])

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = P7P5, fill = P7P5), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(5)) +
  scale_x_continuous(limits = c(0,3)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Hash UMI") +
  ylab("Density") +
  scale_fill_manual(name = "P7 + P5",
                     breaks = c("B02", "E06", "G07"),
                     values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2,3)],
                     labels = c("B + 02",
                                "E + 06",
                                "G + 07"))
ggsave("QC/hash_density_100umi.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(top_to_second_best_ratio), y = timepoint, fill = timepoint), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(2)) +
  scale_x_continuous(limits = c(0,2.5)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Top to Second Best Ratio") +
  ylab("Density") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)])
ggsave("QC/hash_ratio_density_100umi.png", dpi = 600, height = 2, width = 3)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100) %>% arrange(desc(hash_umis)) %>% dplyr::mutate(rank = row_number()), 
       aes(x = log10(rank), y = log10(hash_umis))) +
  geom_point(size = 0.5, stroke = 0, aes(color = P7P5)) +
  geom_hline(linewidth = 0.2, yintercept = log10(10)) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "P7 + P5",
                    breaks = c("B02", "E06", "G07"),
                    values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2,3)],
                    labels = c("B + 02",
                               "E + 06",
                               "G + 07")) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5)) +
  xlab("Log10 CellRank") +
  ylab("Log10 UMI") +
  guides(guides(color = guide_legend(override.aes = list(size=2)))) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom", 
        legend.key.width = unit(0.5,"line"),
        legend.key.height = unit(0.2,"line"),
        text = element_text(size = 6))
ggsave("QC/hash_kneeplot.png", dpi = 600, height = 1.5, width = 2)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100 & hash_umis >= 5),
       aes(x = SHP2_mut, y = percent_mito)) +
  geom_violin(aes(fill = timepoint)) +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100 & hash_umis >= 3 & top_to_second_best_ratio >= 2),
       aes(x = SHP2_mut, y = log10_umi)) +
  geom_violin() +
  monocle3:::monocle_theme_opts() +
  stat_summary(fun.y=median, geom="point", color="black", stroke = 0, size = 1) +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

test <- colData(cds) %>% as.data.frame() %>%
  dplyr::group_by(SHP2_mut) %>%
  dplyr::summarise(median_UMI = median(n_umi))


# ==============================================================================================
cds <- readRDS("SHP2mut_EGFstim_Q510_cds_unfiltered.RDS")

cds <- cds[,colData(cds)$n_umi >= 300]

# predicting and filtering doublets, print matrix file and run scrublet
Matrix::writeMM(t(exprs(cds)), file = "QC/UMI_count_filt.matrix")

scrublet_scores <- read_tsv("QC/scrublet/doublet_scores_EGFRi.txt",
                            col_names = F)

# colData(cds)$doublet_score <- scrublet_scores$X1
# ggplot(colData(cds) %>% as.data.frame, aes(x=doublet_score)) +
#   geom_histogram(color = 'lightpink2', fill = 'lightpink2', binwidth = 0.02) +
#   geom_vline(aes(xintercept=0.50), color = 'black') +
#   annotate("text", label = "Filter = 0.50", x = .62, y = 9000) +
#   monocle3:::monocle_theme_opts() +
#   xlab("Doublet Score") +
#   ylab("Count")
# ggsave('QC/scrublet/doublet_histogram.png',
#        height = 4, width = 6)

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 300 &
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 5 & 
             colData(cds)$top_to_second_best_ratio >= 2]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, "Q510_SHP2_EGFstim_cds_QCfiltered.RDS")
cds <- readRDS("Q510_SHP2_EGFstim_cds_QCfiltered.RDS")

# ================================================================================================================

# Final cell counts
n_umi_per_PCR_tot <- colData(cds) %>% as.data.frame() %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  mutate(P7_row = "Total")

n_umi_per_PCR <- colData(cds) %>% as.data.frame() %>%
  separate(P7, into = c(NA, "P7_row"), sep = 2) %>%
  group_by(P7_row) %>%
  summarise(count = n(), mean_umi_per_cell = mean(n_umi), med_umi_per_cell = median(n_umi), total_umi = sum(n_umi)) %>%
  bind_rows(n_umi_per_PCR_tot)

test <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(count = n(), .groups = "drop") %>%
  # group_by(cell_line) %>%
  summarise(average = mean(count), 
            median = median(count),
            lowest = min(count),
            highest = max(count)) 

test_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(count = n(), med_umi = median(n_umi),.groups = "drop")

ggplot(data = test_plot %>% mutate(count = case_when(count >= 1500 ~ 1500, TRUE ~ count)), 
       aes(x = SHP2_mut, y = count)) +
  facet_wrap(~timepoint) +
  geom_point(aes(color = factor(EGF_dose, levels = c("0", "12.5", "25", "50", "100", "250", "500", "1000"))), 
             position = "jitter") +
  # geom_hline(yintercept = 204, linetype = 2) +
  viridis::scale_color_viridis(name = "EGF (ng/mL)", discrete = T) +
  monocle3:::monocle_theme_opts() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.3,"cm"),
        legend.direction = "horizontal",
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Cell Count Per\nDrug/Dose Condition") +
  guides(color = guide_legend(override.aes = list(size = 2)))
ggsave("QC/cell_count_per_condition_postfilt.png", dpi = 600, height = 2.5, width = 3)

test <- colData(cds) %>% 
  as.data.frame() %>% 
  select(timepoint, hash_plate, hash_well) %>%
  group_by(timepoint, hash_plate, hash_well) %>%
  summarise(count = n()) %>%
  separate(hash_well, into = c("hash_row", "hash_col"), sep = 1, remove = F) %>%
  mutate(hash_col = factor(hash_col, levels = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")),
         hash_row = factor(hash_row, levels = c("H", "G", "F", "E", "D", "C", "B", "A")))

for (cell in c("24hr", "96hr")) {
  ggplot(test %>% filter(timepoint == cell), aes(x = hash_col, y = hash_row, fill = count)) +
    facet_wrap(~hash_plate, scales = "free") +
    geom_tile() +
    viridis::scale_fill_viridis(na.value = "grey30") +
    monocle3:::monocle_theme_opts() +
    scale_x_discrete(position = "top") +
    # scale_y_discrete(expand = c(0,0)) +
    theme(axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(0.5, "line"),
          legend.key.width = unit(0.5, "line"),
          strip.placement = "outside", 
          text = element_text(size = 6))
  ggsave(paste("QC/hash_plate_layout_",cell,".png", sep = ""), height = 2, width = 5, dpi = 600)
}


#========================================================================================
# Correlating replicates for QC
#========================================================================================

expressed_genes <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                            dim(cds)[2]*0.05 ,])

# by line
for (time in c("24hr", "96hr")) {
  for (mut in unique(colData(cds)$SHP2_mut)) {
    cds_temp <- cds[expressed_genes,
                    colData(cds)$timepoint == time & 
                      colData(cds)$SHP2_mut == mut &
                      colData(cds)$hash_plate == "01"]
    rep1 <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
      rownames_to_column(var = "gene")
    cds_temp <- cds[expressed_genes,
                    colData(cds)$timepoint == time & 
                      colData(cds)$SHP2_mut == mut &
                      colData(cds)$hash_plate == "02"]
    rep2 <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp))) %>%
      rownames_to_column(var = "gene")
    
    corr_df <- left_join(rep1, rep2, by = c("gene" = "gene")) %>%
      column_to_rownames(var = "gene")
    corr <- cor(x = corr_df, method = "pearson")
    
    temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
      geom_point(size = 0.2) +
      geom_smooth(method = "lm", 
                  linewidth = 0.3) +
      annotate(geom = "text",
               x = (min(corr_df$rep_1) + 3 * sd(corr_df$rep_1)),
               y = max(corr_df$rep_2),
               label = paste("r =", round(corr[1,2], 3)),
               size = 2) +
      xlab("Replicate 1") +
      ylab("Replicate 2") +
      theme(text = element_text(size = 6)) +
      monocle3:::monocle_theme_opts()
    ggsave(plot = temp, filename = paste("QC/replicate_correlation/normalized_counts_replicate_corr_", time, "_", mut, ".png", sep = ""),
           dpi = 600, height = 2, width = 2)
    
  }
}

for (time in c("24hr", "96hr")) {
  rep1 <- data.frame()
  rep2 <- data.frame()
  for (mut in unique(colData(cds)$SHP2_mut)) {
    cds_temp_1 <- cds[expressed_genes,
                    colData(cds)$timepoint == time & 
                      colData(cds)$SHP2_mut == mut &
                      colData(cds)$hash_plate == "01"]
    rep1_temp <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp_1))) %>%
      rownames_to_column(var = "gene") %>%
      mutate(SHP2_mut = mut)
    rep1 <- bind_rows(rep1, rep1_temp)
    cds_temp_2 <- cds[expressed_genes,
                    colData(cds)$timepoint == time & 
                      colData(cds)$SHP2_mut == mut &
                      colData(cds)$hash_plate == "02"]
    rep2_temp <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp_2))) %>%
      rownames_to_column(var = "gene")%>%
      mutate(SHP2_mut = mut)
    rep2 <- bind_rows(rep2, rep2_temp)
    
  }
  corr_df <- left_join(rep1, rep2, by = c("gene" = "gene",
                                          "SHP2_mut" = "SHP2_mut")) %>%
    unite(col = "gene", sep = "_", gene, SHP2_mut) %>%
    column_to_rownames(var = "gene")
  corr <- cor(x = corr_df, method = "pearson")
  
  temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
    geom_point(size = 0.2) +
    geom_smooth(method = "lm", 
                linewidth = 0.3) +
    annotate(geom = "text",
             x = (min(corr_df$rep_1) + 3 * sd(corr_df$rep_1)),
             y = max(corr_df$rep_2),
             label = paste("r =", round(corr[1,2], 3)),
             size = 2) +
    xlab("Replicate 1") +
    ylab("Replicate 2") +
    theme(text = element_text(size = 6)) +
    monocle3:::monocle_theme_opts()
  ggsave(plot = temp, filename = paste("QC/replicate_correlation/all_genotypes/normalized_counts_replicate_corr_", time,".png", sep = ""),
         dpi = 600, height = 2, width = 2)
  
}

# both timepoints in 1 plot
rep1 <- data.frame()
rep2 <- data.frame()
for (time in c("24hr", "96hr")) {
  for (mut in unique(colData(cds)$SHP2_mut)) {
    cds_temp_1 <- cds[expressed_genes,
                      colData(cds)$timepoint == time & 
                        colData(cds)$SHP2_mut == mut &
                        colData(cds)$hash_plate == "01"]
    rep1_temp <- data.frame("rep_1" = Matrix::rowMeans(normalized_counts(cds_temp_1))) %>%
      rownames_to_column(var = "gene") %>%
      mutate(SHP2_mut = mut) %>%
      mutate(timepoint = time)
    rep1 <- bind_rows(rep1, rep1_temp)
    cds_temp_2 <- cds[expressed_genes,
                      colData(cds)$timepoint == time & 
                        colData(cds)$SHP2_mut == mut &
                        colData(cds)$hash_plate == "02"]
    rep2_temp <- data.frame("rep_2" = Matrix::rowMeans(normalized_counts(cds_temp_2))) %>%
      rownames_to_column(var = "gene")%>%
      mutate(SHP2_mut = mut) %>%
      mutate(timepoint = time)
    rep2 <- bind_rows(rep2, rep2_temp)
    
  }
}
corr_df <- left_join(rep1, rep2, by = c("gene" = "gene",
                                        "SHP2_mut" = "SHP2_mut",
                                        "timepoint" = "timepoint")) %>%
  unite(col = "gene", sep = "_", gene, SHP2_mut, timepoint) %>%
  column_to_rownames(var = "gene")
corr <- cor(x = corr_df, method = "pearson")

temp <- ggplot(data = corr_df, aes(x = rep_1, y = rep_2)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", 
              linewidth = 0.3) +
  annotate(geom = "text",
           x = (min(corr_df$rep_1) + 6 * sd(corr_df$rep_1)),
           y = max(corr_df$rep_2),
           label = paste("r =", round(corr[1,2], 3)),
           size = 3) +
  xlab("Replicate 1") +
  ylab("Replicate 2") +
  theme(text = element_text(size = 8)) +
  monocle3:::monocle_theme_opts()
ggsave(plot = temp, filename = paste("QC/replicate_correlation/all_genotypes/normalized_counts_replicate_corr_alltime.png", sep = ""),
       dpi = 600, height = 2, width = 2)

write_csv(x = corr_df %>% tibble::rownames_to_column(var = "id"), 
          file = "QC/replicate_correlation/all_genotypes/SHP2_Q510_screen_mean_normalized_counts_by_replicate.csv")

# ==========================================================================================
# Getting cell counts, median umi, median hash umi, unique genes per cell line
# ==========================================================================================

cell_count_metrics_tot <- colData(cds) %>% as.data.frame() %>%
  summarise(cell_count = n(), 
            median_umi = median(n_umi),
            median_hash_umi = median(hash_umis),
            mean_express_genes = mean(num_genes_expressed)) %>%
  mutate(cell_line = "Total") %>%
  select(cell_line, everything())

cell_count_metrics <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint) %>%
  summarise(cell_count = n(), 
            median_umi = median(n_umi),
            median_hash_umi = median(hash_umis),
            mean_express_genes = mean(num_genes_expressed))

cell_count_metrics_per_cell_drug_dose <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  summarise(mean_mut_stim_count = mean(cell_count),
            max = max(cell_count),
            min = min(cell_count))

cell_count_metrics_per_drug_dose_tot <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  summarise(mean_mut_stim_count = mean(cell_count),
            max = max(cell_count),
            min = min(cell_count)) %>%
  mutate(cell_line = "Total") %>%
  select(cell_line, everything())

cell_count_metrics_per_drug_dose_tot <- colData(cds) %>% as.data.frame() %>%
  group_by(SHP2_mut) %>%
  summarise(cell_count = n(), mean_umi = mean(n_umi), .groups = "drop")




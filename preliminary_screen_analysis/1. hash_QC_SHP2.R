library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

# adding hash and sci-RNA-seq3 RT barcode sample information
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
  mutate(proportion = hash_umis/sum(hash_umis),total_hash_umis_per_cell = sum(hash_umis)) %>%
  arrange(desc(proportion)) %>%
  mutate(rank = row_number()) %>%
  dplyr::filter(rank %in% c(1,2)) %>%
  mutate(top_to_second_best_ratio = ifelse(sum(rank) > 1, proportion[rank == 1]/proportion[rank == 2], 10)) %>%
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
  geom_hline(linewidth = 0.2, yintercept = log10(350)) +
  monocle3:::monocle_theme_opts() +
  # scale_color_manual(name = "Timepoint",
  #                    breaks = c("24hr", "96hr"),
  #                    values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)]) +
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
ggsave("QC/Kneeplot_allcells.png", dpi = 600, height = 1.25, width = 1.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(hash_umis >= 5)) +
  ggridges::geom_density_ridges(aes(x = log10(n_umi), y = timepoint, fill = timepoint), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(300)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 UMI") +
  ylab("Density") +
  scale_x_continuous(breaks = seq(0,5,1)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)])
ggsave("QC/umi_density.png", dpi = 600, height = 2, width = 2.5)

saveRDS(cds, "SHP2mut_EGFstim_cds_unfiltered.RDS")

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(hash_umis), y = timepoint, fill = timepoint), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(10)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Hash UMI") +
  ylab("Density") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)])
ggsave("QC/hash_density_100umi.png", dpi = 600, height = 2, width = 2.5)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100)) +
  ggridges::geom_density_ridges(aes(x = log10(top_to_second_best_ratio), y = timepoint, fill = timepoint), 
                                alpha = 0.5, 
                                rel_min_height = 0.001,
                                position = position_nudge(y = -0.5),
                                show.legend = FALSE) +
  geom_vline(xintercept = log10(3)) +
  monocle3:::monocle_theme_opts() +
  xlab("Log10 Top to Second Best Ratio") +
  ylab("Density") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)])
ggsave("QC/hash_ratio_density_100umi.png", dpi = 600, height = 2, width = 3)

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100) %>% arrange(desc(hash_umis)) %>% dplyr::mutate(rank = row_number()), 
       aes(x = log10(rank), y = log10(hash_umis))) +
  geom_point(size = 0.5, stroke = 0, aes(color = timepoint)) +
  geom_hline(linewidth = 0.2, yintercept = log10(10)) +
  monocle3:::monocle_theme_opts() +
  scale_color_manual(name = "Timepoint",
                     breaks = c("24hr", "96hr"),
                     values = RColorBrewer::brewer.pal(n = 4, "Set1")[c(1,2)]) +
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

ggplot(colData(cds) %>% as.data.frame() %>% filter(n_umi >= 100 & hash_umis >= 5 & top_to_second_best_ratio >= 2.5),
       aes(x = SHP2_mut, y = log10_umi)) +
  geom_violin() +
  monocle3:::monocle_theme_opts() +
  viridis::scale_fill_viridis(option = "B", discrete = TRUE) +
  theme(axis.text.x = element_blank())

# ==============================================================================================
cds <- readRDS("SHP2mut_EGFstim_cds_unfiltered.RDS")

# Final QC cutoffs
cds <- cds[,colData(cds)$n_umi >= 300 &
             !is.na(colData(cds)$hash_umis) &
             colData(cds)$hash_umis >= 5 & 
             colData(cds)$top_to_second_best_ratio >= 2.5]

cds <- estimate_size_factors(cds)
cds <- detect_genes(cds)

source("cell_cycle.R")
cc.genes <- readRDS("cc.genes.RDS")
cds <- estimate_cell_cycle(cds, g1s_markers = cc.genes$s.genes, g2m_markers = cc.genes$g2m.genes)

saveRDS(cds, "preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

# ================================================================================================================
cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

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
  group_by(SHP2_mut) %>%
  summarise(count = n(), .groups = "drop") %>%
  summarise(average = mean(count), 
            median = median(count),
            lowest = min(count),
            highest = max(count)) 

test <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(count = n(), .groups = "drop") %>%
  summarise(average = mean(count), 
            median = median(count),
            lowest = min(count),
            highest = max(count)) 

test_plot <- colData(cds) %>% as.data.frame() %>%
  group_by(timepoint, SHP2_mut, EGF_dose) %>%
  summarise(count = n(), med_umi = median(n_umi),.groups = "drop")

#========================================================================================
# Correlating replicates for QC
#========================================================================================

expressed_genes <- row.names(rowData(cds)[Matrix::rowSums(exprs(cds) > 0) > 
                                            dim(cds)[2]*0.05 ,])

# each individual time and mutant
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

# combine mutants, keeping time separate
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


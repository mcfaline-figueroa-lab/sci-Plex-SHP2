library(devtools)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(monocle3)
library(ComplexHeatmap)
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

colData(cds)$timepoint <- factor(colData(cds)$timepoint, levels = c("24hr", "96hr"))

colData(cds)$EGF_log10_dose <- log10(as.double(colData(cds)$EGF_dose) + 0.1)

colData(cds)$replicate <- case_when(
  colData(cds)$hash_plate == "01" ~ 1,
  colData(cds)$hash_plate == "02" ~ 2
)

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.01,]$id

# filter SHP2_mut KO and Q510E (transfection control) to better identify significant changes for other mut
cds <- cds[,!colData(cds)$SHP2_mut %in% c("KO", "Q510E")]

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

saveRDS(treatment_diff_test, "SHP2_toWT_full_treatment_diff_test.RDS")
treatment_diff_test <- readRDS("SHP2_toWT_full_treatment_diff_test.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(q_value <= 0.05)

treatment_diff_test_SHP2_plot <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  mutate(plot_color = case_when(
    abs(normalized_effect) >= 0.1 & q_value <= 0.01 ~ "red",
    TRUE ~ "black"
  )) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("WT", "KO",
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
           plot_color == "red" & gene_short_name == "PTPN11" ~ gene_short_name,
           TRUE ~ NA
         )),
       aes(x = normalized_effect, y = plot_q_value)) +
  facet_grid(timepoint~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3,
                           min.segment.length = 0.01, 
                           force_pull = 0.7,
                           segment.size = 0.2
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
ggsave("SHP2_volcano_plot_highlight.png",
       dpi = 600, height = 2.5, width = 6.5)


for (select_gene in c("EGFR", "KRAS", "NRAS", "MAGED1", "PTPN11", "SH3BP2", "VGF", "JUNB")) {
  ggplot(treatment_diff_test_SHP2_plot %>% 
           mutate(plot_gene = case_when(
             gene_short_name == select_gene ~ gene_short_name,
             TRUE ~ NA
           )) %>%
           mutate(plot_color = case_when(
             !is.na(plot_gene) ~ "red",
             TRUE ~ "black"
           )) %>%
           arrange(plot_color) %>%
           mutate(plot_color = factor(plot_color)),
         aes(x = normalized_effect, y = plot_q_value)) +
    facet_grid(timepoint~SHP2_mut) +
    geom_point(aes(color = plot_color),
               size = 0.2,
               show.legend = F) +
    ggrepel::geom_text_repel(aes(label = plot_gene), 
                             size = 1.3,
                             max.overlaps = 8,
                             min.segment.length = 0.1, 
                             force_pull = 0,
                             segment.size = 0.1
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
  ggsave(paste0("SHP2_volcano_plot_",select_gene,"_highlight.png"),
         dpi = 600, height = 2.5, width = 6.5)
}

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
  facet_grid(timepoint~SHP2_mut) +
  geom_point(aes(color = plot_color),
             size = 0.2,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3,
                           min.segment.length = 0.1, 
                           force_pull = 0.7,
                           segment.size = 0.1
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
       dpi = 600, height = 2.5, width = 6.5)

ggplot(treatment_diff_test_SHP2_plot %>% 
         mutate(plot_gene = case_when(
           gene_short_name %in% pathway_genes ~ gene_short_name,
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
             size = 0.3,alpha = 0.4,
             stroke = 0.01,
             show.legend = F) +
  ggrepel::geom_text_repel(aes(label = plot_gene), 
                           size = 1.3,
                           min.segment.length = 0.05, 
                           force_pull = 0.4,
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
ggsave("SHP2_volcano_plot_RTK_highlight.png",
       dpi = 600, height = 2.5, width = 3.5)


treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  group_by(timepoint, SHP2_mut, direction) %>%
  dplyr::summarise(count = n()) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K",
                                                "Q510E", "Q510K")))

ggplot(treatment_diff_test_SHP2 %>% mutate(count = case_when(
  count >= 160 ~ 160,
  TRUE ~ count
)), 
       aes(x = SHP2_mut, y = count)) +
  facet_wrap(~timepoint) +
  geom_bar(stat = "identity", 
           aes(fill = direction),
           position = "dodge") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(n = 3, "Set1")[c(2,1)]) +
  scale_y_continuous(breaks = c(0,40,80,120,160), labels = c("0", "40", "80", "120", "160+")) +
  xlab("SHP2 Mutant") +
  ylab("Count") +
  theme(text = element_text(size = 6),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(linewidth = 0.1),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(0,0,l = -0.5,0)) +
  guides(fill = guide_legend(title = "FDR < 5%\nabs(beta coeff) > 0.1")) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_mutant_significant_genes.png",
       dpi = 600, height = 1.5, width = 3)

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  filter(!SHP2_mut %in% c("KO", "Q510E")) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  distinct(gene_short_name, direction) %>%
  group_by(direction) %>%
  dplyr::summarise(total_count = n())

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  filter(!SHP2_mut %in% c("KO", "Q510E")) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  group_by(SHP2_mut, direction) %>%
  dplyr::summarise(count = n(), .groups = "drop") %>%
  group_by(direction) %>%
  summarise(median_count = median(count))

# ==============================================================================
# Plotting individual genes
# ==============================================================================

colData(cds)$EGF_dose_group <- case_when(
  colData(cds)$EGF_dose == 0 ~ "No EGF",
  colData(cds)$EGF_dose %in% c(12.5, 25, 50) ~ "Low",
  colData(cds)$EGF_dose %in% c(100, 250, 500, 1000) ~ "High"
)

test <- cds[rowData(cds)$gene_short_name %in% c("EGFR"),
           colData(cds)$timepoint == "24hr"]
colData(test)$condition <- paste0(colData(test)$SHP2_mut, "_", colData(test)$EGF_dose_group)

plot_percent_cells_positive(cds = cds[rowData(cds)$gene_short_name %in% c("EGFR"),],
                            group_cells_by = "SHP2_mut") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ==============================================================================
# Analyzing intersections of DEGs between SHP2 mutants
# ==============================================================================

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  ))

upset.list <- list()
for (mut in unique(treatment_diff_test_SHP2$SHP2_mut) %>% sort()) {
  upset.list[[mut]] <- treatment_diff_test_SHP2 %>% 
    filter(SHP2_mut == mut & direction == "Upregulated") %>% 
    pull(gene_short_name)
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) > 3], pt_size = unit(1, "mm"), lwd = 0.6,
                                        comb_order = order(comb_size(intersection[comb_size(intersection) > 3])), 
                                        set_order = c("KO",
                                                      "T42A", "T52S",
                                                      "E76K", "R138Q",
                                                      "E139D", "Y279C",
                                                      "T468M", "T507K",
                                                      "Q510E", "Q510K"),
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

png("SHP2_mut_DEG_upreg_intersection.png",width=4.5,height=2,units="in",res=1200)
draw(intersect_plot)
dev.off()

intersection <- intersection[comb_size(intersection) >= 4]
gene_sets <- sapply(comb_name(intersection), function(nm) extract_comb(intersection, nm))
test <- set_name(intersection)
test_names <- data.frame(x = names(gene_sets)) %>%
  separate(x, into = as.character(c(1:11)), sep = c(1,2,3,4,5,6,7,8,9,10,11)) %>%
  mutate(`1` = ifelse(`1` == 1, test[1],""),
         `2` = ifelse(`2` == 1, test[2],""),
         `3` = ifelse(`3` == 1, test[3],""),
         `4` = ifelse(`4` == 1, test[4],""),
         `5` = ifelse(`5` == 1, test[5],""),
         `6` = ifelse(`6` == 1, test[6],""),
         `7` = ifelse(`7` == 1, test[7],""),
         `8` = ifelse(`8` == 1, test[8],""),
         `9` = ifelse(`9` == 1, test[9],""),
         `10` = ifelse(`10` == 1, test[10],""),
         `11` = ifelse(`11` == 1, test[11],"")) %>%
  unite(col = "name", sep = " ", 1:11)

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

# saveRDS(gene_sets_final, "SHP2_mut_DEGup_intersection.RDS")

upset.list <- list()
for (mut in unique(treatment_diff_test_SHP2$SHP2_mut) %>% sort()) {
  upset.list[[mut]] <- treatment_diff_test_SHP2 %>% 
    filter(SHP2_mut == mut & direction == "Downregulated") %>% 
    pull(id)
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

intersect_plot <- ComplexHeatmap::UpSet(intersection[comb_size(intersection) > 3], pt_size = unit(1, "mm"), lwd = 0.6,
                                        comb_order = order(comb_size(intersection[comb_size(intersection) > 3])), 
                                        set_order = c("KO",
                                                      "T42A", "T52S",
                                                      "E76K", "R138Q",
                                                      "E139D", "Y279C",
                                                      "T468M", "T507K",
                                                      "Q510E", "Q510K"),
                                        top_annotation = upset_top_annotation(intersection[comb_size(intersection) > 3],
                                                                              # ylim = c(0, 500),
                                                                              add_numbers = FALSE,
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
                                        column_title = "SHP2mut DEG Intersection (downregulated, FDR < 0.01)",
                                        column_title_gp = gpar(fontsize = 8, fontface = "bold"))

png("SHP2_mut_DEG_downreg_intersection.png",width=4.5,height=2,units="in",res=1200)
draw(intersect_plot)
dev.off()

intersection <- intersection[comb_size(intersection) >= 4]
gene_sets <- sapply(comb_name(intersection), function(nm) extract_comb(intersection, nm))
test <- set_name(intersection)
test_names <- data.frame(x = names(gene_sets)) %>%
  separate(x, into = as.character(c(1:11)), sep = c(1,2,3,4,5,6,7,8,9,10,11)) %>%
  mutate(`1` = ifelse(`1` == 1, test[1],""),
         `2` = ifelse(`2` == 1, test[2],""),
         `3` = ifelse(`3` == 1, test[3],""),
         `4` = ifelse(`4` == 1, test[4],""),
         `5` = ifelse(`5` == 1, test[5],""),
         `6` = ifelse(`6` == 1, test[6],""),
         `7` = ifelse(`7` == 1, test[7],""),
         `8` = ifelse(`8` == 1, test[8],""),
         `9` = ifelse(`9` == 1, test[9],""),
         `10` = ifelse(`10` == 1, test[10],""),
         `11` = ifelse(`11` == 1, test[11],"")) %>%
  unite(col = "name", sep = " ", 1:11)

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

# Grab DEG test information for each element in list
gene_sets_final <- list()
for (new_element in names(gene_sets)) {
  new_element_split <- stringr::str_split_1(string = new_element, pattern = "_")
  gene_set_df_temp <- treatment_diff_test_SHP2 %>%
    filter(id %in% gene_sets[[new_element]]) %>%
    filter(SHP2_mut %in% new_element_split) %>%
    select(-id, -status, -estimate, -std_err, -test_val, -p_value, -model_component) %>%
    as.data.frame()
  gene_sets_final[[new_element]] <- gene_set_df_temp
}

# saveRDS(gene_sets_final, "SHP2_mut_DEGdown_intersection.RDS")

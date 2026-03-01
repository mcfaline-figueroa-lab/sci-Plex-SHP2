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


expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.05,]$id

treatment_diff_test <- fit_models(cds = cds[expressed_genes,], 
                                  model_formula_str = "~ timepoint + SHP2_mut + EGF_log10_dose + replicate",
                                  expression_family = "quasipoisson",
                                  cores = 1, 
                                  verbose = TRUE)
treatment_diff_test <- coefficient_table(treatment_diff_test) %>% dplyr::select(-model,-model_summary)
treatment_diff_test$SHP2_mut <- stringr::str_remove(treatment_diff_test$term,
                                                    "SHP2_mut")

saveRDS(treatment_diff_test, "SHP2_full_treatment_diff_test.RDS")
treatment_diff_test <- readRDS("SHP2_full_treatment_diff_test.RDS")

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
ggsave("SHP2_volcano_plot_highlight.png",
       dpi = 600, height = 2.5, width = 3.5)

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
       dpi = 600, height = 2.5, width = 3.5)

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
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  group_by(SHP2_mut, direction) %>%
  dplyr::summarise(count = n()) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("KO", "WT",
                                                "T42A", "T52S",
                                                "E76K", "R138Q",
                                                "E139D", "Y279C",
                                                "T468M", "T507K",
                                                "Q510E", "Q510K")))

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
  guides(fill = guide_legend(title = "FDR < 1%\nabs(beta coeff) > 0.1")) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_mutant_significant_genes.png",
       dpi = 600, height = 1.5, width = 3)

# ==============================================================================
# Analyzing intersections of DEGs between SHP2 mutants
# ==============================================================================

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.1 & q_value <= 0.01) %>%
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
                                        set_order = c("WT",
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

saveRDS(gene_sets_final, "SHP2_mut_DEGup_intersection.RDS")

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
                                        set_order = c("WT",
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
  new_element_split <- str_split_1(string = new_element, pattern = "_")
  gene_set_df_temp <- treatment_diff_test_SHP2 %>%
    filter(id %in% gene_sets[[new_element]]) %>%
    filter(SHP2_mut %in% new_element_split) %>%
    select(-id, -status, -estimate, -std_err, -test_val, -p_value, -model_component) %>%
    as.data.frame()
  gene_sets_final[[new_element]] <- gene_set_df_temp
}

saveRDS(gene_sets_final, "SHP2_mut_DEGdown_intersection.RDS")

# ================================================================================
# Overlap of DEGs for mock or WT vs DMSO
# ================================================================================

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) > 0.25 & q_value < 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  )) %>%
  filter(SHP2_mut %in% c("WT", "Q510E")) %>%
  filter(direction == "Upregulated") %>%
  filter(timepoint == "24hr") %>%
  mutate(SHP2_mut = case_when(
    SHP2_mut == "Q510E" ~ "Mock Transfected",
    TRUE ~ SHP2_mut
  ))

upset.list <- list()
for (mut in unique(treatment_diff_test_SHP2$SHP2_mut) %>% sort()) {
  upset.list[[mut]] <- treatment_diff_test_SHP2 %>% 
    filter(SHP2_mut == mut & direction == "Upregulated") %>% 
    pull(gene_short_name)
}

intersection <- ComplexHeatmap::make_comb_mat(upset.list,
                                              mode = "distinct")

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
names(gene_sets) <- c("Mock Transfected_WT", "Mock Transfected", "WT")

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

test <- gene_sets_final[["Mock Transfected_WT"]] %>%
  select(gene_short_name, normalized_effect, SHP2_mut, timepoint, direction) %>%
  pivot_wider(id_cols = c(gene_short_name,
                          timepoint,
                          direction),
              names_from = SHP2_mut,
              values_from = normalized_effect)

readr::write_csv(x = test,
                 file = "SHP2_WT_mock_upregulated_intersection.csv",
                 col_names = T)

gene_universe <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.25 & q_value <= 0.05) %>%
  select(gene_short_name) %>%
  distinct(gene_short_name) %>%
  pull()

source("GSA_helper_functions.r")
source("loadGSCSafe.R")

hallmarksGSC <- loadGSCSafe(file = "h.all.v6.0.symbols.gmt",
                            type = "gmt")

hallmarksGSC_length <- data.frame("gene_set" = names(hallmarksGSC[["gsc"]])) %>%
  tibble::rownames_to_column("row")

GSAhyper_df_final <- data.frame()
contingency_tab_list <- list()
for (temp_name in c("WT", "Mock Transfected_WT", "Mock Transfected")) {
  gene_test <- gene_sets_final[[temp_name]] %>%
    select(gene_short_name) %>% 
    distinct(gene_short_name) %>% 
    pull()
  
  test <- piano::runGSAhyper(genes = gene_test, 
                             gsc = hallmarksGSC, 
                             adjMethod = "BH", 
                             universe = gene_universe)
  
  GSAhyper_df <- as.data.frame(test$pvalues) %>%
    tibble::rownames_to_column()
  colnames(GSAhyper_df) <- c("gene_set","p_value")
  
  for (i in GSAhyper_df$gene_set) {
    
    hallmarkGSC_index <- hallmarksGSC_length %>%
      mutate(row = as.double(row)) %>%
      filter(gene_set == i) %>%
      pull(row)
    
    contingency_tab_list[[paste0(temp_name)]][[i]] <- test$contingencyTable[[c(hallmarkGSC_index)]]
  }
  
  GSAhyper_df_final <- bind_rows(GSAhyper_df_final, 
                                 GSAhyper_df %>%
                                   mutate(intersection = temp_name))
  
}

GSAhyper_df_final$q_value <- p.adjust(GSAhyper_df_final$p_value)

ggplot(GSAhyper_df_final %>%
         filter(q_value <= 0.01) %>%
         mutate(intersection = case_when(
           intersection == "WT" ~ "WT Alone",
           intersection == "Mock Transfected" ~ "Mock Trasnfected\nAlone",
           TRUE ~ "WT /\nMock Transfected"
         )) %>%
         mutate(gene_set = stringr::str_remove(gene_set,
                                               "HALLMARK_")),
       aes(y = gene_set, x = intersection)) +
  geom_point(aes(color = -log10(q_value),
                 size = -log10(q_value),
                 alpha = intersection)) +
  viridis::scale_color_viridis(option = "D") +
  scale_alpha_manual(values = c(1,1,0),
                     breaks = c("WT Alone", "WT /\nMock Transfected", "Mock Trasnfected\nAlone")) +
  scale_size_continuous(range = c(1,4)) +
  xlab("DEG Intersection") +
  ylab("Hallmark Gene Set") +
  guides(size = "none",
         alpha = "none") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.margin = margin(t = -5, l = -100)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_WT_mock_upregulated_intersection_GSEA.png",
       dpi = 900, height = 3, width = 3)

readr::write_csv(x = GSAhyper_df_final %>%
                   filter(q_value <= 0.01) %>%
                   mutate(intersection = case_when(
                     intersection == "WT" ~ "WT Alone",
                     TRUE ~ "WT /\nMock Transfected"
                   )) %>%
                   mutate(log10_FDR = -log10(q_value)),
                 file = "SHP2_WT_mock_upregulated_intersection_GSEA_res_for_printing.csv",
                 col_names = T)


# ==============================================================================
# Aggregate expression of union SHP2 mutant genes
# ==============================================================================
treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(abs(normalized_effect) >= 0.05 & q_value <= 0.05) %>%
  mutate(direction = case_when(
    normalized_effect >= 0 ~ "Upregulated",
    TRUE ~ "Downregulated"
  ))

DEG_union <- treatment_diff_test_SHP2 %>%
  distinct(id) %>%
  pull()

colData(cds)$condition <- paste0(colData(cds)$timepoint,"_",
                                 colData(cds)$SHP2_mut)

cds_agg <- cds[DEG_union,
               colData(cds)$EGF_dose %in% c(100, 250, 500, 1000)]

cell_group_df <- tibble::tibble(cell=colData(cds_agg)$Cell, 
                                cell_group=colData(cds_agg)$condition)

agg_mat <- aggregate_gene_expression(cds_agg, 
                                     gene_group_df = NULL,
                                     cell_group_df = cell_group_df, 
                                     norm_method = "log",
                                     pseudocount = 1, 
                                     scale_agg_values = T) %>%
  as.data.frame() %>%
  select(contains("24hr"))

colnames(agg_mat) <- gsub("24hr_","",colnames(agg_mat))
colnames(agg_mat) <- gsub("Q510E","Mock Transfected",colnames(agg_mat))

agg_mat[agg_mat > 2] <- 2
agg_mat[agg_mat < -2] <- -2

hm <- ComplexHeatmap::Heatmap(matrix = agg_mat %>% as.matrix(), 
                              name = "Z-Scored\nExpression",
                              heatmap_legend_param = list(title = "Z-Scored\nExpression",
                                                          title_gp = gpar(fontface = "bold", 
                                                                          fontsize = 8),
                                                          labels_gp = gpar(fontsize = 7),
                                                          grid_width = unit(0.3, "cm"),
                                                          legend_height = unit(1.8, "cm")),
                              show_column_names = TRUE, 
                              column_names_gp = gpar(fontsize = 8), 
                              column_names_rot = 45,
                              column_split = 6,
                              row_split = 8,
                              column_title = NULL,
                              cluster_columns = TRUE, 
                              cluster_rows = TRUE,
                              row_title = "Gene", 
                              row_title_gp = gpar(fontsize = 6),
                              show_row_names = F,
                              row_dend_width = unit(5,"mm"),
                              column_dend_height = unit(5, "mm"),
                              use_raster = F)
draw(hm)

hm <- ComplexHeatmap::Heatmap(matrix = agg_mat %>% as.matrix() %>% t(), 
                              name = "Z-Scored\nExpression",
                              heatmap_legend_param = list(title = "Z-Scored\nExpression",
                                                          title_gp = gpar(fontface = "bold", 
                                                                          fontsize = 8),
                                                          labels_gp = gpar(fontsize = 7),
                                                          grid_height = unit(0.3, "cm"),
                                                          legend_width = unit(1.8, "cm"), 
                                                          direction = "horizontal",
                                                          title_position = "leftcenter"),
                              show_row_names = TRUE, 
                              row_names_gp = gpar(fontsize = 8), 
                              row_split = 6,
                              column_title = "Module %s",
                              column_title_side = "bottom",
                              column_title_rot = 45,
                              column_split = 10,
                              row_gap = unit(1.5, "mm"),
                              cluster_columns = TRUE, 
                              cluster_rows = TRUE,
                              row_title = NULL, 
                              column_title_gp = gpar(fontsize = 8),
                              show_column_names = F,
                              row_dend_width = unit(5,"mm"),
                              column_dend_height = unit(5, "mm"),
                              use_raster = F)
draw(hm,
     heatmap_legend_side = "bottom")

pdf("SHP2_24hr_EGFhigh_aggregate_expression_DEG_union.pdf",width=6.5,height=3.5,compress = F)
ht <- draw(hm,
           heatmap_legend_side = "bottom"); col_clusters <- column_order(ht); row_clusters <- row_order(ht)
dev.off()

col_order_heatmap_mat <- data.frame("column" = colnames(agg_mat %>% as.matrix() %>% t())) %>%
  rownames_to_column(var = "row_numb")
col_order_heatmap_mat_final <- data.frame()
for(name in 1:length(col_clusters)) {
  col_order_heatmap_mat_temp <- col_order_heatmap_mat %>%
    filter(row_numb %in% col_clusters[[name]]) %>%
    mutate(`Gene Cluster` = name)
  col_order_heatmap_mat_final <- bind_rows(col_order_heatmap_mat_final, 
                                           col_order_heatmap_mat_temp)
}

col_order_heatmap_mat_final_named <- col_order_heatmap_mat_final %>%
  select(-row_numb) %>%
  dplyr::rename(id = column) %>%
  left_join(rowData(cds) %>% as.data.frame() %>%
              select(-num_cells_expressed),
            by = c("id"))

write_csv(x = col_order_heatmap_mat_final_named %>%
            mutate(Module = `Gene Cluster`), 
          file = "SHP2_24hr_EGFhigh_aggregate_expression_DEG_union_gene_modules.csv")

gene_module_cluster_summary <- col_order_heatmap_mat_final_named %>%
  group_by(`Gene Cluster`) %>%
  summarise(count = n())

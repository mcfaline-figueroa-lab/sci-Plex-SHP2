library(tidyverse)
library(monocle3)
library(ComplexHeatmap)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

treatment_diff_test <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

treatment_diff_test$q_value <- p.adjust(treatment_diff_test$p_value, method = "BH")

treatment_diff_test_SHP2 <- treatment_diff_test %>%
  filter(grepl(pattern = "SHP2_mut", term)) %>%
  filter(q_value <= 0.05) %>%
  filter(abs(normalized_effect) >= 0.05)

transfection_control_genes <- treatment_diff_test_SHP2 %>%
  filter(SHP2_mut == "Q510E")

core_cutoff_df <- data.frame()
for (time in c("24hr", "96hr")) {
  transfection_control_genes <- treatment_diff_test_SHP2 %>%
    filter(SHP2_mut == "Q510E") %>%
    filter(timepoint == time)
  for (i in 1:10) {
    temp_up <- treatment_diff_test_SHP2 %>%
      filter(SHP2_mut != "Q510E") %>%
      filter(timepoint == time) %>%
      filter(normalized_effect > 0) %>%
      filter(!id %in% (transfection_control_genes %>% 
               filter(normalized_effect > 0) %>%
               pull(id))) %>%
      group_by(id, gene_short_name) %>%
      summarise(count = n(), .groups = "drop") %>%
      filter(count >= i) %>%
      summarise(total_count = n()) %>%
      mutate(timepoint = time,
             cutoff = i,
             direction = "Upregulated")
    
    temp_down <- treatment_diff_test_SHP2 %>%
      filter(SHP2_mut != "Q510E") %>%
      filter(timepoint == time) %>%
      filter(normalized_effect < 0) %>%
      filter(!id %in% (transfection_control_genes %>% 
               filter(normalized_effect < 0) %>%
               pull(id))) %>%
      group_by(id, gene_short_name) %>%
      summarise(count = n(), .groups = "drop") %>%
      filter(count >= i) %>%
      summarise(total_count = n()) %>%
      mutate(timepoint = time,
             cutoff = i,
             direction = "Downregulated")
    
    core_cutoff_df <- bind_rows(core_cutoff_df, temp_up, temp_down)
    
  }
  
}

ggplot(core_cutoff_df %>%
         mutate(total_count = ifelse(total_count >= 1000, 1000, total_count)) %>%
         mutate(plot_count = ifelse(total_count >= 1000, NA, as.character(total_count))), 
       aes(x = cutoff, y = total_count)) +
  facet_grid(timepoint ~ direction) +
  geom_text(aes(label = plot_count),
            size = 1,
            nudge_y = 75,
            nudge_x = 0.4) +
  geom_point(size = 0.5) +
  geom_line(linewidth = 0.3) +
  geom_hline(data = hline_6mut,
             aes(yintercept = total_count),
             linetype = 2,
             linewidth = 0.2) +
  scale_x_continuous(breaks = seq(1,10,1),
                     limits = c(1,11)) +
  xlab("SHP2 variant # cutoff") +
  scale_y_continuous(limits = c(0,1100),
                     breaks = seq(0,1000,200),
                     labels = c("0", "200", "400", "600", "800", ">1000")) +
  ylab("# core SHP2-dependent\ntranscripts") +
  ggtitle("Core Transcriptome Size (FDR < 5%)") +
  theme(plot.title = element_text(size = 7,
                                  hjust = 0.5),
        text = element_text(size = 6),
        strip.text = element_text(size = 6),
        axis.ticks = element_line(linewidth = unit(0.25, "mm")),
        axis.ticks.length = unit(0.75, "mm")) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_core_interactome_cutoff_counts.png",
       dpi = 600, height = 2, width = 2.75)


write_csv(core_cutoff_df, file = "SHP2_core_transcriptome_cutoff_counts.csv")

core_cutoff_df <- data.frame()
for (time in c("24hr", "96hr")) {
  transfection_control_genes <- treatment_diff_test_SHP2 %>%
    filter(SHP2_mut == "Q510E") %>%
    filter(timepoint == time)
  for (i in 5) {
    temp_up <- treatment_diff_test_SHP2 %>%
      filter(SHP2_mut != "Q510E") %>%
      filter(timepoint == time) %>%
      filter(normalized_effect > 0) %>%
      filter(!id %in% (transfection_control_genes %>% 
                         filter(normalized_effect > 0) %>%
                         pull(id))) %>%
      group_by(id, gene_short_name) %>%
      summarise(count = n(), .groups = "drop") %>%
      filter(count >= i) %>%
      mutate(timepoint = time,
             cutoff = i,
             direction = "Upregulated")
    
    temp_down <- treatment_diff_test_SHP2 %>%
      filter(SHP2_mut != "Q510E") %>%
      filter(timepoint == time) %>%
      filter(normalized_effect < 0) %>%
      filter(!id %in% (transfection_control_genes %>% 
                         filter(normalized_effect < 0) %>%
                         pull(id))) %>%
      group_by(id, gene_short_name) %>%
      summarise(count = n(), 
                .groups = "drop") %>%
      filter(count >= i) %>%
      mutate(timepoint = time,
             cutoff = i,
             direction = "Downregulated")
    
    core_cutoff_df <- bind_rows(core_cutoff_df, temp_up, temp_down)
    
  }
  
}

core_cutoff_df_for_printing <- core_cutoff_df %>%
  select(-cutoff) %>%
  dplyr::rename(SHP2_variant_DEG_count = count)

write_csv(core_cutoff_df_for_printing,
          "SHP2_core_transcriptome_by_time_DEGdirection.csv")

test_summary <- core_cutoff_df_for_printing %>%
  group_by(timepoint, direction) %>%
  summarise(final_count = n())

# =================================================================================
# Core transcriptome characterization
# =================================================================================

core_SHP2_transcriptome <- read_csv("SHP2_core_transcriptome_by_time_DEGdirection.csv") %>%
  filter(SHP2_variant_DEG_count == 5)

core_24hr_up <- core_SHP2_transcriptome %>%
  filter(timepoint == "24hr" & direction == "Upregulated")

core_96hr_up <- core_SHP2_transcriptome %>%
  filter(timepoint == "96hr" & direction == "Upregulated")

core_24hr_dn <- core_SHP2_transcriptome %>%
  filter(timepoint == "24hr" & direction == "Downregulated")

core_96hr_dn <- core_SHP2_transcriptome %>%
  filter(timepoint == "96hr" & direction == "Downregulated")

# Select gene universe (tested genes)
gene_universe <- treatment_diff_test %>%
  select(gene_short_name) %>%
  distinct()


# download helper functions from sci-Plex github and gmt files from online source
source("GSA_helper_functions.r")
source("loadGSCSafe.R")

hallmarksGSC <- loadGSCSafe(file = "h.all.v6.0.symbols.gmt",
                            type = "gmt")
oncogenicGSC <- loadGSCSafe(file = "c6.all.v2023.2.Hs.symbols.gmt",
                            type = "gmt")
kinaseupGSC <- loadGSCSafe(file = "gene_set_library_up_crisp.gmt",
                           type = "gmt")
kinasedownGSC <- loadGSCSafe(file = "gene_set_library_dn_crisp.gmt",
                             type = "gmt")


GSAhyper_df_final <- data.frame()
for (time in c("24hr", "96hr")) {
  for (dir in c("Upregulated", "Downregulated")) {
    test_df <- core_SHP2_transcriptome %>%
      filter(timepoint == time & direction == dir)
    test <- piano::runGSAhyper(genes = test_df$gene_short_name, 
                               gsc = hallmarksGSC, 
                               adjMethod = "BH", 
                               universe = gene_universe$gene_short_name)
    
    GSAhyper_df <- as.data.frame(test$p.adj) %>%
      rownames_to_column()
    colnames(GSAhyper_df) <- c("gene_set","q_value")
    GSAhyper_df <- GSAhyper_df %>%
      mutate(core_SHP2 = paste0("core_SHP2_",time,"_",dir),
             gsc_set = "HALLMARKS") %>%
      filter(q_value <= 0.25)
    
    GSAhyper_df_final <- bind_rows(GSAhyper_df_final, GSAhyper_df)
    
  }
  
}

write_csv(file = "core_SHP2_transcriptome_hypergeometric_GSA_results.csv",
          x = GSAhyper_df_final, col_names = T)


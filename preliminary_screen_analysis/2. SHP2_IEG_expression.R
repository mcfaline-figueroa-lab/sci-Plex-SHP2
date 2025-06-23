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
source("calculate_aggreg_expression.R")

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

colData(cds)$EGF_dose_group <- case_when(
  colData(cds)$EGF_dose == 0 ~ "No EGF",
  colData(cds)$EGF_dose %in% c(12.5, 25, 50) ~ "Low",
  colData(cds)$EGF_dose %in% c(100, 250, 500, 1000) ~ "High"
)

colData(cds)$EGF_dose_group <- factor(colData(cds)$EGF_dose_group,
                                      levels = c("No EGF", "Low", "High"))

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.01,]$id

# filter SHP2_mut KO and transfection control to better identify significant changes for other mut
cds <- cds[,colData(cds)$SHP2_mut %in% c("WT", "KO")]

# Genes from https://pubmed.ncbi.nlm.nih.gov/17575275/
IEG_genes <- read_csv("Immediate_Early_Genes.csv") %>%
  dplyr::rename(Gene = 1) %>%
  mutate(Gene = toupper(Gene))

colData(cds)$IEG <- calculate_aggregate_expression_score(cds,
                                                         signature_genes = IEG_genes$Gene)

colData(cds)$PTPN11 <- calculate_aggregate_expression_score(cds,
                                                            signature_genes = c("PTPN11"))

ggplot(colData(cds) %>% as.data.frame(),
       aes(x = PTPN11)) +
  geom_density(outline.type = "full") +
  monocle3:::monocle_theme_opts()

ggplot(colData(cds) %>% as.data.frame() %>%
         mutate(PTPN11_bin = 
                  case_when(
                    PTPN11 <= 0  ~ "Not Expressing",
                    TRUE ~ "Expressing"
                  )) %>%
         mutate(PTPN11_bin = factor(PTPN11_bin, 
                                    levels = c("Not Expressing", "Expressing"))) %>%
         mutate(SHP2_mut_PTPN11_bin = case_when(
           SHP2_mut == "WT" ~ PTPN11_bin,
           TRUE ~ "KO"
           )) %>%
         mutate(SHP2_mut_PTPN11_bin = factor(SHP2_mut_PTPN11_bin, levels = c("KO", "Not Expressing", "Expressing"))),
       aes(x = SHP2_mut_PTPN11_bin, y = IEG)) +
  facet_wrap(~EGF_dose_group) +
  # facet_grid(timepoint~EGF_dose_group) +
  geom_violin(aes(fill = SHP2_mut_PTPN11_bin),
              show.legend = F,
              scale = "area",
              linewidth = 0.15) +
  stat_summary(geom = "point",
               fun = "median",
               size = 0.75) + 
  # scale_fill_manual(values = c("lightpink", "grey60")) +
  xlab("PTPN11") +
  ylab("Intermediate Early Response\nGene Expression") +
  ggtitle("Expression of IEG program\nwithin SHP2 WT") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8, hjust = 0.5)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_WT_IEG_expression_by_PTPN11_bin_KO_incl.png",
       dpi = 600, height = 2.25, width = 2.25)

WT_IEG_data_for_printing <- colData(cds) %>% as.data.frame() %>%
  mutate(PTPN11_bin = 
           case_when(
             PTPN11 <= 0  ~ "Not Expressing",
             TRUE ~ "Expressing"
           )) %>%
  mutate(PTPN11_bin = factor(PTPN11_bin, 
                             levels = c("Not Expressing", "Expressing"))) %>%
  mutate(SHP2_mut_PTPN11_bin = case_when(
    SHP2_mut == "WT" ~ PTPN11_bin,
    TRUE ~ "KO"
  )) %>%
  mutate(SHP2_mut_PTPN11_bin = factor(SHP2_mut_PTPN11_bin, levels = c("KO", "Not Expressing", "Expressing"))) %>%
  select(Cell, timepoint, SHP2_mut, EGF_dose, EGF_dose_group, IEG_expression = IEG, PTPN11_expression = PTPN11, PTPN11_bin, x_axis_for_plot = SHP2_mut_PTPN11_bin)


# reload cds object to look at IEG program expression of mutants
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

colData(cds)$EGF_dose_group <- case_when(
  colData(cds)$EGF_dose == 0 ~ "No EGF",
  colData(cds)$EGF_dose %in% c(12.5, 25, 50) ~ "Low",
  colData(cds)$EGF_dose %in% c(100, 250, 500, 1000) ~ "High"
)

colData(cds)$EGF_dose_group <- factor(colData(cds)$EGF_dose_group,
                                      levels = c("No EGF", "Low", "High"))

colData(cds)$IEG <- calculate_aggregate_expression_score(cds,
                                                         signature_genes = IEG_genes$Gene)

colData(cds)$PTPN11 <- calculate_aggregate_expression_score(cds,
                                                            signature_genes = c("PTPN11"))


ggplot(colData(cds) %>% as.data.frame() %>%
         mutate(PTPN11_bin = 
                  case_when(
                    PTPN11 <= 0  ~ "Not Expressing",
                    TRUE ~ "Expressing"
                  )) %>%
         mutate(PTPN11_bin = factor(PTPN11_bin, 
                                    levels = c("Not Expressing", "Expressing"))),
       aes(x = PTPN11_bin, y = IEG)) +
  facet_grid(EGF_dose_group ~ SHP2_mut) +
  geom_violin(aes(fill = PTPN11_bin),
              show.legend = F,
              scale = "area",
              linewidth = 0.15) +
  stat_summary(geom = "point",
               fun = "median",
               size = 0.75) + 
  scale_fill_manual(values = c("lightpink", "grey60")) +
  xlab("PTPN11") +
  ylab("Intermediate Early Response\nGene Expression") +
  ggtitle("Expression of IEG program\nwithin SHP2 WT") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8, hjust = 0.5)) +
  monocle3:::monocle_theme_opts()


ggplot(colData(cds) %>% as.data.frame() %>%
         # filter(timepoint == "24hr") %>%
         mutate(PTPN11_bin = 
                  case_when(
                    PTPN11 <= 0  ~ "Not Expressing",
                    TRUE ~ "Expressing"
                  )) %>%
         mutate(PTPN11_bin = factor(PTPN11_bin, 
                                    levels = c("Not Expressing", "Expressing"))) %>%
         filter(PTPN11_bin == "Expressing") %>%
         filter(!SHP2_mut %in% c("KO", "Q510E")),
       aes(x = SHP2_mut, y = IEG)) +
  facet_wrap(~EGF_dose_group) +
  geom_violin(aes(fill = EGF_dose_group),
              show.legend = F,
              scale = "area",
              linewidth = 0.15) +
  stat_summary(geom = "point",
               fun = "median",
               size = 0.75) + 
  # scale_fill_manual(values = c("lightpink", "grey60")) +
  xlab("PTPN11") +
  ylab("Intermediate Early Response\nGene Expression") +
  ggtitle("Expression of IEG program\nonly SHP2-expressing cells") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8, hjust = 0.5)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_mut_IEG_expression_only_PTPN11_expressing.png",
       dpi = 600, height = 2.25, width = 4)

ggplot(colData(cds) %>% as.data.frame() %>%
         mutate(PTPN11_bin = 
                  case_when(
                    PTPN11 <= 0  ~ "Not Expressing",
                    TRUE ~ "Expressing"
                  )) %>%
         mutate(PTPN11_bin = factor(PTPN11_bin, 
                                    levels = c("Not Expressing", "Expressing"))) %>%
         filter(PTPN11_bin == "Expressing") %>%
         filter(!SHP2_mut %in% c("KO", "Q510E")),
       aes(x = EGF_dose_group, y = IEG)) +
  facet_wrap(~SHP2_mut, nrow = 2) +
  geom_violin(aes(fill = EGF_dose_group),
              show.legend = F,
              scale = "area",
              linewidth = 0.15) +
  stat_summary(geom = "point",
               fun = "median",
               size = 0.75) + 
  # scale_fill_manual(values = c("lightpink", "grey60")) +
  xlab("EGF Exposure Group") +
  ylab("Intermediate Early Response\nGene Expression") +
  ggtitle("Expression of IEG program\nwithin SHP2 expressing cells") +
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 8, hjust = 0.5)) +
  monocle3:::monocle_theme_opts()
ggsave("SHP2_mut_IEG_expression_only_PTPN11_expressing_by_mut.png",
       dpi = 600, height = 3, width = 2.75)

all_mut_IEG_data_for_printing <- colData(cds) %>% as.data.frame() %>%
  # filter(timepoint == "24hr") %>%
  mutate(PTPN11_bin = 
           case_when(
             PTPN11 <= 0  ~ "Not Expressing",
             TRUE ~ "Expressing"
           )) %>%
  mutate(PTPN11_bin = factor(PTPN11_bin, 
                             levels = c("Not Expressing", "Expressing"))) %>%
  select(Cell, timepoint, SHP2_mut, EGF_dose, EGF_dose_group, IEG_expression = IEG, PTPN11_expression = PTPN11, PTPN11_bin)



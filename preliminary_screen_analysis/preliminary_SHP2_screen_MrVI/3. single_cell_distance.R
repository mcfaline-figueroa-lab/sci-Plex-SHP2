library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

SHP2_mut_factor <- c("KO", "WT",
                     "T42A", "T52S",
                     "E76K", "R138Q",
                     "E139D", "Y279C",
                     "T468M", "T507K",
                     "Q510E", "Q510K")

temp_df <- read_csv("SHP2_24hr_cell_dists_namedosesample_attention_singlecell.csv") %>%
  filter(target_drug %in% c("24hr_WT_0.0")) %>%
  separate(sample_drug, into = c("sample_timepoint", "sample_SHP2_mut", "sample_EGF_dose"),
           convert = T, remove = F, sep = "_") %>%
  separate(target_drug, into = c("target_timepoint", "target_SHP2_mut", "target_EGF_dose"),
           convert = T, remove = F, sep = "_")

temp_df_plot <- temp_df %>%
  filter(sample_SHP2_mut %in% c("WT", "R138Q", "E76K", "T42A")) %>%
  dplyr::rename(EGF_dose = sample_EGF_dose,
                SHP2_mut = sample_SHP2_mut) %>%
  mutate(EGF_dose = case_when(
    EGF_dose == 12.5 ~ round(EGF_dose, 1),
    TRUE ~ round(EGF_dose, 0)
  )) %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    EGF_dose %in% c("100", "250", "500", "1000") ~ "High"
  )) %>%
  mutate(EGF_log_dose = case_when(
    EGF_dose != 0 ~ round(log10(as.double(EGF_dose)), 2),
    TRUE ~ 0
    )) %>%
  mutate(EGF_log_dose = factor(EGF_log_dose)) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = c("WT", "R138Q", "E76K", "T42A"))) %>%
  mutate(EGF_dose_group = factor(EGF_dose_group, 
                                 levels = c("No EGF", "Low", "High")))

temp_df_plot_mean <- temp_df_plot %>%
  dplyr::group_by(SHP2_mut, EGF_dose_group) %>%
  dplyr::summarise(mean_distance = mean(`Counterfactual Distance`),
                   median_distance = median(`Counterfactual Distance`))

temp_df_plot_hline <- temp_df_plot_mean %>%
  filter(EGF_dose_group == "No EGF")

ggplot(data = temp_df_plot, 
       aes(x = EGF_dose_group, y = `Counterfactual Distance`, fill = EGF_dose_group)) +
  facet_wrap(~SHP2_mut, nrow = 1) +
  # geom_boxplot(outlier.shape = NA) +
  # scale_y_continuous(limits = c(0, 1.25)) +
  geom_violin(show.legend = F, linewidth = 0.25) +
  geom_point(data = temp_df_plot_mean,
             aes(y = mean_distance),
             show.legend = F, size = 0.5) +
  ylab("Counterfactual Distance\nto WT 0ng/ml EGF") +
  # geom_hline(data = temp_df_plot_hline,
  #            aes(yintercept = mean_distance)) +
  geom_hline(yintercept = temp_df_plot_hline %>% filter(SHP2_mut == "R138Q") %>% pull(mean_distance),
             linewidth = 0.25, linetype = 2) +
  theme(text = element_text(size = 6),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()
ggsave("single_cell_distances_to_WTNoEGF.png",
       dpi = 600, height = 1.5, width = 2.5)

temp_df_for_csv <- temp_df %>%
  dplyr::rename(EGF_dose = sample_EGF_dose,
                SHP2_mut = sample_SHP2_mut) %>%
  mutate(EGF_dose = case_when(
    EGF_dose == 12.5 ~ round(EGF_dose, 1),
    TRUE ~ round(EGF_dose, 0)
  )) %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    EGF_dose %in% c("100", "250", "500", "1000") ~ "High"
  )) %>%
  rownames_to_column(var = "cell") %>%
  select(cell, sample_name = sample_drug, 
         SHP2_mut, EGF_dose, 
         EGF_dose_group, target_sample = target_drug,
         `Counterfactual Distance`)

write_csv(x = temp_df_for_csv,
          "SHP2_24hr_MrVI_single_cell_dists_to_WTnoEGF.csv")  

# ===========================================================================
# Counterfactual distances to a variant's own No EGF sample
# ===========================================================================

temp_df <- read_csv("SHP2_24hr_cell_dists_namedosesample_attention_singlecell.csv") %>%
  separate(sample_drug, into = c("sample_timepoint", "sample_SHP2_mut", "sample_EGF_dose"),
           convert = T, remove = F, sep = "_") %>%
  separate(target_drug, into = c("target_timepoint", "target_SHP2_mut", "target_EGF_dose"),
           convert = T, remove = F, sep = "_") %>%
  filter(sample_SHP2_mut == target_SHP2_mut & target_EGF_dose == 0)

temp_df_plot <- temp_df %>%
  dplyr::rename(EGF_dose = sample_EGF_dose,
                SHP2_mut = sample_SHP2_mut) %>%
  mutate(EGF_dose = case_when(
    EGF_dose == 12.5 ~ round(EGF_dose, 1),
    TRUE ~ round(EGF_dose, 0)
  )) %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    EGF_dose %in% c("100", "250", "500", "1000") ~ "High"
  )) %>%
  mutate(EGF_log_dose = case_when(
    EGF_dose != 0 ~ round(log10(as.double(EGF_dose)), 2),
    TRUE ~ 0
  )) %>%
  mutate(EGF_log_dose = factor(EGF_log_dose)) %>%
  mutate(SHP2_mut = factor(SHP2_mut, levels = SHP2_mut_factor)) %>%
  mutate(EGF_dose_group = factor(EGF_dose_group, 
                                 levels = c("No EGF", "Low", "High")))

temp_df_plot_mean <- temp_df_plot %>%
  dplyr::group_by(SHP2_mut, EGF_dose_group) %>%
  dplyr::summarise(mean_distance = mean(`Counterfactual Distance`),
                   median_distance = median(`Counterfactual Distance`))

temp_df_plot_hline <- temp_df_plot_mean %>%
  filter(EGF_dose_group == "No EGF")

ggplot(data = temp_df_plot %>%
         filter(EGF_dose_group != "No EGF") %>%
         filter(!SHP2_mut %in% c("KO", "Q510E")), 
       aes(x = EGF_dose_group, y = `Counterfactual Distance`, fill = EGF_dose_group)) +
  facet_wrap(~SHP2_mut, nrow = 1) +
  scale_y_continuous(limits = c(0.1, 1.5)) +
  geom_violin(show.legend = F, linewidth = 0.25) +
  geom_point(data = temp_df_plot_mean %>%
               filter(EGF_dose_group != "No EGF") %>%
               filter(!SHP2_mut %in% c("KO", "Q510E")),
             aes(y = median_distance),
             show.legend = F, size = 0.5) +
  ylab("Counterfactual Distance\nto 0ng/ml EGF") +
  theme(text = element_text(size = 8),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()

ggplot(data = temp_df_plot %>%
         filter(EGF_dose_group != "No EGF") %>%
         filter(!SHP2_mut %in% c("KO", "Q510E")) %>%
         mutate(EGF_dose = log10(EGF_dose)), 
       aes(x = EGF_dose, y = `Counterfactual Distance`)) +
  facet_wrap(~SHP2_mut, nrow = 1) +
  geom_smooth(method = "loess") +
  # geom_smooth(method = "glm", formula = y ~ log(x)) +
  ylab("Counterfactual Distance\nto 0ng/ml EGF") +
  theme(text = element_text(size = 8),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()

ggplot(data = temp_df_plot %>%
         filter(EGF_dose %in% c(12.5, 25, 50, 100, 500, 1000)) %>%
         filter(!SHP2_mut %in% c("KO", "Q510E")), 
       aes(x = as.factor(EGF_dose), y = `Counterfactual Distance`, fill = EGF_dose_group)) +
  facet_wrap(~SHP2_mut, nrow = 1) +
  # geom_boxplot(outlier.shape = NA) +
  # scale_y_continuous(limits = c(0, 1.25)) +
  geom_violin(show.legend = F, linewidth = 0.25) +
  # geom_point(data = temp_df_plot_mean,
  #            aes(y = median_distance),
  #            show.legend = F, size = 0.5) +
  ylab("Counterfactual Distance\nto 0ng/ml EGF") +
  # geom_hline(data = temp_df_plot_hline,
  #            aes(yintercept = mean_distance)) +
  # geom_hline(yintercept = temp_df_plot_hline %>% filter(SHP2_mut == "R138Q") %>% pull(mean_distance),
  #            linewidth = 0.25, linetype = 2) +
  theme(text = element_text(size = 6),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()
  
ggplot(data = temp_df_plot %>%
         filter(EGF_dose != 0), 
       aes(x = as.factor(EGF_dose), y = `Counterfactual Distance`,
           group = SHP2_mut)) +
  facet_wrap(~SHP2_mut, nrow = 1) +
  geom_smooth(method = "loess") +
  # geom_point(data = temp_df_plot_mean,
  #            aes(y = median_distance),
  #            show.legend = F, size = 0.5) +
  ylab("Counterfactual Distance\nto 0ng/ml EGF") +
  # geom_hline(data = temp_df_plot_hline,
  #            aes(yintercept = mean_distance)) +
  # geom_hline(yintercept = temp_df_plot_hline %>% filter(SHP2_mut == "R138Q") %>% pull(mean_distance),
  #            linewidth = 0.25, linetype = 2) +
  theme(text = element_text(size = 6),
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  monocle3:::monocle_theme_opts()

temp_df_for_csv <- temp_df %>%
  dplyr::rename(EGF_dose = sample_EGF_dose,
                SHP2_mut = sample_SHP2_mut) %>%
  mutate(EGF_dose = case_when(
    EGF_dose == 12.5 ~ round(EGF_dose, 1),
    TRUE ~ round(EGF_dose, 0)
  )) %>%
  mutate(EGF_dose_group = case_when(
    EGF_dose == "0" ~ "No EGF",
    EGF_dose %in% c("12.5", "25", "50") ~ "Low",
    EGF_dose %in% c("100", "250", "500", "1000") ~ "High"
  )) %>%
  filter(!SHP2_mut %in% c("KO", "Q510E")) %>%
  rownames_to_column(var = "cell") %>%
  select(cell, sample_name = sample_drug, 
         SHP2_mut, EGF_dose, 
         EGF_dose_group, target_sample = target_drug,
         `Counterfactual Distance`)

write_csv(x = temp_df_for_csv,
          "SHP2_24hr_MrVI_single_cell_dists_to_respective_noEGF.csv")  


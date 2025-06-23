library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("Q510_SHP2_EGFstim_cds_QCfiltered.RDS")

cds <- detect_genes(cds)
expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.05,] %>%
  as.data.frame() %>%
  pull(id)

cds_subset <- cds[expressed_genes,
                  !colData(cds)$SHP2_mut %in% c("HEK293", "R4A-R5A", "KO")]

# write to mtx and csv files
Matrix::writeMM(t(exprs(cds_subset)), 
                file = "SHP2_count_all.mtx")

write_csv(x = colData(cds_subset) %>% as.data.frame(),"SHP2_colData.csv")
write_csv(x = rowData(cds_subset) %>% as.data.frame(),"SHP2_rowData.csv")






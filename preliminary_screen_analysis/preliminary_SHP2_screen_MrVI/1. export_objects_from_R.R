library(tidyverse)
library(monocle3)
this.dir <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(this.dir))
rm(this.dir)

cds <- readRDS("preliminary_SHP2_EGFstim_cds_QCfiltered.RDS")

expressed_genes <-  rowData(cds)[rowData(cds)$num_cells_expressed >= nrow(cds)*0.01,] %>% 
  as.data.frame()

cds_subset <- cds[expressed_genes$id,]

# predicting and filtering doublets, print matrix file and run scrublet
Matrix::writeMM(t(exprs(cds_subset)), 
                file = "SHP2_count_all.mtx")

write_csv(x = colData(cds_subset) %>% as.data.frame(),"SHP2_colData.csv")
write_csv(x = rowData(cds_subset) %>% as.data.frame(),"SHP2_rowData.csv")






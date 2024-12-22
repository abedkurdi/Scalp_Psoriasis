library(Seurat) 
library(metap) 
library(MAST) 
library(scDblFinder) 
library(ggplot2) 
library(cowplot) 
library(patchwork) 
library(stringr) 
library(RColorBrewer) 
library(ComplexHeatmap) 
library(circlize) 
library(tidyverse) 
library(concaveman) 
library(ggforce) 
library(monocle3) 
library(dplyr) 
library(Rmisc) 
library(ggrepel) 
library(raster) 
library(vegan)
library(tools)
library("writexl")
library(openxlsx)
library(harmony)

# load in the data
merged_samples <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_sct.RDS")

# set the resolution 0.1
Idents(merged_samples) <- "SCT_snn_res.0.1"

# annotate the cells at resolution 0.1
merged_samples$cluster <- as.character(merged_samples$SCT_snn_res.0.1)

merged_samples$cluster[merged_samples$cluster == 0] <- "Lymphocyte"
merged_samples$cluster[merged_samples$cluster == 1] <- "Lymphocyte"
merged_samples$cluster[merged_samples$cluster == 2] <- "Lymphocyte"
merged_samples$cluster[merged_samples$cluster == 3] <- "APC"
merged_samples$cluster[merged_samples$cluster == 4] <- "APC"
merged_samples$cluster[merged_samples$cluster == 5] <- "Mast"
merged_samples$cluster[merged_samples$cluster == 6] <- "Fibroblast"

# checkpoint 1: save the annotated object
saveRDS(merged_samples, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_sct_and_annotated.RDS")

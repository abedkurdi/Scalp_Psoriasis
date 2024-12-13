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
merged_samples <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated.RDS")

# set the resolution 0.1
Idents(merged_samples) <- "RNA_snn_res.0.1"

# annotate the cells at resolution 0.1
merged_samples$cluster <- as.character(merged_samples$RNA_snn_res.0.1)
merged_samples$supercluster <- as.character(merged_samples$RNA_snn_res.0.1)

merged_samples$cluster[merged_samples$cluster == 0] <- "T"
merged_samples$cluster[merged_samples$cluster == 1] <- "Treg"
merged_samples$cluster[merged_samples$cluster == 2] <- "CTL"
merged_samples$cluster[merged_samples$cluster == 3] <- "DC"
merged_samples$cluster[merged_samples$cluster == 4] <- "Mac"
merged_samples$cluster[merged_samples$cluster == 5] <- "Mast"
merged_samples$cluster[merged_samples$cluster == 6] <- "Mac/Tet"
merged_samples$cluster[merged_samples$cluster == 7] <- "Fibroblast"
merged_samples$cluster[merged_samples$cluster == 8] <- "Neutrophils/Mac"
merged_samples$cluster[merged_samples$cluster == 9] <- "B"
merged_samples$cluster[merged_samples$cluster == 10] <- "cycling"
merged_samples$cluster[merged_samples$cluster == 11] <- "NK"
merged_samples$cluster[merged_samples$cluster == 12] <- "Mono"

merged_samples$supercluster[merged_samples$supercluster == 0] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 1] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 2] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 3] <- "APC"
merged_samples$supercluster[merged_samples$supercluster == 4] <- "APC"
merged_samples$supercluster[merged_samples$supercluster == 5] <- "Mast"
merged_samples$supercluster[merged_samples$supercluster == 6] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 7] <- "Fibroblast"
merged_samples$supercluster[merged_samples$supercluster == 8] <- "APC"
merged_samples$supercluster[merged_samples$supercluster == 9] <- "B"
merged_samples$supercluster[merged_samples$supercluster == 10] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 11] <- "Lymphocyte"
merged_samples$supercluster[merged_samples$supercluster == 12] <- "APC"

# checkpoint 1: save the annotated object
saveRDS(merged_samples, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_and_annotated.RDS")

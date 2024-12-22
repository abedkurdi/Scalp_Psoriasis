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

output_dir <- "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_sct/"
samples_dir <- "/special_projects/55_scalp_psoriasis_analysis_jeff/samples/"
metadata_df <- read.csv("/special_projects/55_scalp_psoriasis_analysis_jeff/metadata.txt", header=TRUE, sep='\t', stringsAsFactors=FALSE)
stress <- read.table("/special_projects/55_scalp_psoriasis_analysis_jeff/coregene_df-FALSE-v3.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)

setwd(samples_dir)
all_samples <- list()
all_samples_data <- list()
all_samples_sce <- list()

# read in the samples, make the Seurat objects, remove doublets and plot featurescatter figures
for(f in dir(samples_dir)){
	print(f)

	# read in the samples and remove doublets
	if(length(dir(paste0(samples_dir,"/",f,"/")))==3){
	all_samples_data[[f]] <- Read10X(paste0(samples_dir,"/",f))
	all_samples_data[[f]] <- CreateSeuratObject(counts = all_samples_data[[f]])
	} else {
		all_samples_data[[f]] <- Read10X_h5(paste0(samples_dir,"/",f,"/filtered_feature_bc_matrix.h5"))
		if(!is.null(names(all_samples_data[[f]]))){
			all_samples_data[[f]] <- CreateSeuratObject(all_samples_data[[f]][["Gene Expression"]])
		} else {
			all_samples_data[[f]] <- CreateSeuratObject(counts = all_samples_data[[f]])
		}
	}
	
	# plot the rank of cells according to the number of unique genes
	cat("Plotting the rank of cells per the unique number of genes\n")
	
	pdf(paste0(output_dir,"/",f,"_rank_of_cells_per_unique_genes_before_any_filtering_step.pdf"))
	print(ggplot(data=all_samples_data[[f]]@meta.data, aes(x=reorder(rownames(all_samples_data[[f]]@meta.data),nFeature_RNA), y=nFeature_RNA))+
    geom_point()+ggtitle(f)+geom_hline(yintercept=6000)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()))
	dev.off()
	
	all_samples_sce[[f]] <- as.SingleCellExperiment(all_samples_data[[f]])
	all_samples_sce[[f]] <- scDblFinder(all_samples_sce[[f]])
	print(table(call=all_samples_sce[[f]]$scDblFinder.class))

	all_samples[[f]] <- as.Seurat(all_samples_sce[[f]], project=f)
	DefaultAssay(all_samples[[f]]) <- "RNA"

	all_samples[[f]][["percent.mt"]] <- PercentageFeatureSet(all_samples[[f]], pattern="^MT-", assay="RNA")
	all_samples[[f]][["percent.ribo"]] <- PercentageFeatureSet(all_samples[[f]], pattern="^RPS|^RPL|^MRPS|^MRPL", assay="RNA")
	all_samples[[f]][["percent.stress"]] <- PercentageFeatureSet(all_samples[[f]], features = intersect(stress$gene_symbol, rownames(all_samples[[f]])))
	
	all_samples[[f]] <- subset(all_samples[[f]], scDblFinder.class == "singlet")
	
	# add the metadata
	tmp_df <- subset(metadata_df, sample == f)
	print(paste0("adding metadata for sample ",f))

	for(m in colnames(tmp_df)){
		print(paste0("adding ",m))
		all_samples[[f]][[m]] <- tmp_df[[m]]
	}
}


# according to the upper scatterplots, I set the thresholds
print("applying filters")
for(f in names(all_samples)){
	print(f)
	all_samples[[f]] <- subset(all_samples[[f]], 
	subset = nFeature_RNA > 200 & 
		nFeature_RNA < 6000 & 
		percent.mt < 20 & 
		percent.ribo < 40 & 
		percent.stress < 25)
}

# merge all samples in one object
## rename cells to avoid identical barcodes from different samples
for(s in names(all_samples)){
    print(s)
    all_samples[[s]] <- RenameCells(all_samples[[s]], new.names=paste0(s,"_",colnames(all_samples[[s]])))
}

## Run SCTransform
all_samples <- lapply(X = all_samples, 
                       FUN = SCTransform, 
                       return.only.var.genes = FALSE, ncells = 1000, conserve.memory = TRUE)


## merge all
merged_samples <- merge(all_samples[[1]], y = all_samples[-1], add.cell.ids = NULL)
gc()

# (!) checkpoint One
saveRDS(merged_samples, file="/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_filtered_with_metadata_sct.RDS")

# FindVariableFeatures, set VariableFeatures and RunPCA
var.features <- SelectIntegrationFeatures(object.list = all_samples, nfeatures = 3000)
VariableFeatures(merged_samples) <- var.features
merged_samples <- RunPCA(merged_samples, verbose = FALSE)

# Run Harmony
merged_samples <- merged_samples %>% 
    RunHarmony(c("sample","condition","chemistry","kit_10x"), plot_convergence = TRUE, assay.use='RNA')

pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_sct/elbow_plot.pdf")
ElbowPlot(merged_samples, reduction = "harmony", ndims=30) # PC24 is a good cut-off
dev.off()

## another "quantitative" way to choose PC cut-off
pct <- merged_samples@reductions$harmony@stdev / sum(merged_samples@reductions$harmony@stdev) * 100
cum <- cumsum(pct)
co1 <- which(cum > 90 & pct < 5)[1]
co1
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1),  decreasing = T)[1] + 1 # last point where change of % of variation is more than 0.1%.
co2
pcs <- min(co1, co2)
pcs # 39 is the cut-off according to this method

### According to Bayesian Information Criterion method, <> is the cutoff (skipped)
#selDim <- FindNumberPC(merged_samples, dim2Test = 40, returnBIC = TRUE)
#selDim$selPCA

# (!) checkpoint Two
saveRDS(merged_samples, file="/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_no_clustering_sct.RDS")

dims <- 1:12
neighbors <- 100
metric <- 'cosine'
prune <- 1/15
seed <- 1

merged_samples <- merged_samples %>% 
  FindNeighbors(reduction = 'harmony', dims = dims, k.param = neighbors, annoy.metric = metric, prune.SNN = prune) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 1), random.seed = seed) %>% 
  RunUMAP(reduction = 'harmony', dims = dims, n.neighbors = neighbors, metric = metric, seed.use = seed) 


# generate the marker genes for all resolutions
merged_samples <- PrepSCTFindMarkers(merged_samples)
OUT <- createWorkbook()
resolutions <- unique(colnames(merged_samples@meta.data)[grepl("SCT_snn",colnames(merged_samples@meta.data))],)

for(res in resolutions[c(1,5)]){
	print(res)
	tryCatch({
		addWorksheet(OUT, res)
		Idents(merged_samples) <- res
		markers <- FindAllMarkers(merged_samples, only.pos=TRUE, logfc.threshold=0.25, min.pct=0.25, test.use="MAST")
		writeData(OUT, sheet = res, x = markers)
	}, error=function(e){message(paste0("We have error!",res))})
}

saveWorkbook(OUT, "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_sct/marker_genes_resolutions_0.1_and_0.5.xlsx")

# generate the UMAPs at different resolutions
for(res in resolutions[c(1,5)]){
	Idents(merged_samples) <- res
	pdf(paste0(output_dir,"UMAP_res_",res,".pdf"), width=8, height=8)
	print(DimPlot(merged_samples, label=TRUE, repel=TRUE, raster=FALSE, group.by=res)+NoLegend())
	dev.off()
}



# (!) checkpoint Three
saveRDS(merged_samples, file="/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_sct.RDS")

#-------------------------------------------------------------------------------------------------------------------------------------------------------#


# Dot plots
merged_samples <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_sct.RDS")

rashx_markers <- c("CD3D",
                  "CCR7", "SELL", "KLF2","CD69", #Tn/Tcm/Tmm
                  "ITGAE","CXCR6", #Trm
                  "CD4","FOXP3", "TIGIT", "CTLA4", "IL2RA", #Treg
                  "CD8A","CD8B","GZMB", "NKG7", "CCL5", #CTL
                  "KLRB1", "GATA3","PTGDR2", "PRF1","KLRD1","GNLY", #NK/ILC
                  "TNFRSF4", "CD96","TNFRSF18","PRDM1","PDCD1","LAG3", #ac or ex
                  "BATF", "SNHG12", "ZFAS1", #Tet
                  "TRAT1","RORA","IL7R",  #not specific
                  "HLA-DRA","HLA-DRB1","CD83","IDO1", #APC
                  "CD207","EPCAM", #LC
                  "CD68","CEBPB","FCER1G","C1QB","C1QC","CD163","NR4A1","NR4A2","KLF4","CD14","S100A9", #Mac
                  "MS4A7","LYZ","SERPINA1", "IL23A","IL1B", "CXCL3", #Mono
                  "CD1C", "CLEC10A","CLEC9A", "XCR1","FSCN1", "LAMP3",  #DC
                  "CSF3R", "FCGR3B",  "NCF2",  #neutrophils "CXCR2","CD177", "CEACAM8", "ELANE",
                  "IGKC","JCHAIN","CD79A","MS4A1","CD19", "MS4A2", "IGHG1", "IGHG4", #B
                  "TPSB2","TPSAB1", #Mast
                  "MKI67","TOP2A", "CENPF", "UBE2C", #cycling
                  "THBD","SIRPA","F13A1", #not specific
                  "ITGA4","NCR1","IL17A","IL17F","IL23R", #not specific
                   'PECAM1',  #'COL1A1',  #stromal
                  'KRT1','KRT10','KRT14',#'KRT5', #keratinocyte
                  'COL1A1','COL3A1','COL1A2',   #fibroblast
                  'TAGLN','ACTA2','RGS5',#'MYL9',  #smooth muscle cells
                  'KRT19','KRT7','CLDN10',#'TMPRSS2',  #epithelial
                  'DCT','TYRP1','MLANA',#'PMEL',  #Melanocytes
                 'AQP1','CLDN5','PLVAP','EMCN',  #'SELE', #Vasular endothelial cells
                  'CCL21','TFF3','LYVE1' #'STAB2'  #lymphatic endothelial cells
                  )

## dot plot - resolution 0.1
Idents(merged_samples) <- "SCT_snn_res.0.1"

pdf(paste0(output_dir,"dot_plot_rashx_markers_resolution_0.1.pdf"), width=26, height=6)
DotPlot(merged_samples, features = rashx_markers, cluster.ident=TRUE)+RotatedAxis()
dev.off()

## dot plot - resolution 0.5
Idents(merged_samples) <- "SCT_snn_res.0.5"

pdf(paste0(output_dir,"dot_plot_rashx_markers_resolution_0.5.pdf"), width=26, height=8)
DotPlot(merged_samples, features = rashx_markers, cluster.ident=TRUE)+RotatedAxis()
dev.off()







#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------------------------------------------------------#

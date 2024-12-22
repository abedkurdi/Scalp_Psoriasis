# perform subclustering
merged_samples <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_sct_and_annotated.RDS")
seed <- 123
set.seed(seed)

# variables to use in Harmony
batch <- c('condition', 'sample', 'chemistry', 'kit_10x')

# get some metadata
merged_samples$supercluster <- merged_samples$cluster
celltypes <- merged_samples@meta.data %>% dplyr::select(sample, sample_id, supercluster)

apcs_barcodes <- celltypes[celltypes$supercluster == "APC",] %>% rownames()
lymphocytes_barcodes <- celltypes[celltypes$supercluster == "Lymphocyte",] %>% rownames()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# APCs and Lymphocytes LOCAL subclustering
apcs <- subset(merged_samples, subset = supercluster == "APC")
apcs <- SCTransform(apcs, verbose=FALSE, method="glmGamPoi")
apcs <- RunPCA(apcs)

lymphocyte <- subset(merged_samples, subset = supercluster == "Lymphocyte")
lymphocyte <- SCTransform(lymphocyte, verbose=FALSE, method="glmGamPoi")
lymphocyte <- RunPCA(lymphocyte)

## check the elbow plot
ElbowPlot(apcs)
ElbowPlot(lymphocyte)

## Run Harmony
max_harmony_dims <- 20
apcs <- apcs %>% RunHarmony(group.by.vars = batch,
                            dims.use = 1:max_harmony_dims, ncores=16)

apcs %>% ElbowPlot(ndims = max_harmony_dims, reduction = 'harmony')

lymphocyte <- lymphocyte %>% RunHarmony(group.by.vars = batch,
                            dims.use = 1:max_harmony_dims, ncores=16)

lymphocyte %>% ElbowPlot(ndims = max_harmony_dims, reduction = 'harmony')

## recluster after Harmony
dims <- 1:max_harmony_dims

### define the "recluster" function
recluster <- function(seuratobj, 
                      assay = 'SCT',
                      reduction = 'harmony',
                      dims = 1:15,
                      seed = 1234,
                      n.neighbors = 100,
                      n.epochs = 200,
                      min.dist = 0.3,
                      k.param = n.neighbors,
                      metric = 'cosine',
                      algorithm = 1,
                      resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 1),
                      default_ident = 'SCT_snn_res.0.5') {
  
    ### RunUMAP, FindNeighbors, and FindClusters 
  seuratobj <- seuratobj %>% 
    FindNeighbors(reduction = reduction, 
                  dims = dims, 
                  assay = assay,
                  annoy.metric = metric,
                  k.param = k.param) %>% 
    FindClusters(random.seed = seed,
                 algorithm = algorithm,
                 resolution = resolution) %>% 
    RunUMAP(reduction = reduction, 
            dims = dims, 
            assay = assay, 
            seed.use = seed,
            n.neighbors = n.neighbors,
            n.epochs = n.epochs,
            min.dist = min.dist)
  
  ### Fix cluster factor levels
  clusterings <- colnames(seuratobj@meta.data) %>% str_subset('_res')
  
  for(i in clusterings) {
    clusters <- seuratobj@meta.data[[i]]
    seuratobj[[i]] <- factor(clusters, levels = levels(clusters) %>% as.numeric() %>% sort())
  }
  
  Idents(seuratobj) <- default_ident
  seuratobj$seurat_clusters <- seuratobj[[default_ident]]
  
  seuratobj
}

## do the reclustering  at resolutions 0.1, 0.2, 0.3, 0.4, 0.5 and 1
apcs <- apcs %>% recluster(dims = dims)
lymphocyte <- lymphocyte %>% recluster(dims = dims)

apcs$SCT_snn_res.0.1 <- as.character(apcs$SCT_snn_res.0.1)
apcs$SCT_snn_res.0.2 <- as.character(apcs$SCT_snn_res.0.2)
apcs$SCT_snn_res.0.3 <- as.character(apcs$SCT_snn_res.0.3)
apcs$SCT_snn_res.0.4 <- as.character(apcs$SCT_snn_res.0.4)
apcs$SCT_snn_res.0.5 <- as.character(apcs$SCT_snn_res.0.5)
apcs$SCT_snn_res.1 <- as.character(apcs$SCT_snn_res.1)

lymphocyte$SCT_snn_res.0.1 <- as.character(lymphocyte$SCT_snn_res.0.1)
lymphocyte$SCT_snn_res.0.2 <- as.character(lymphocyte$SCT_snn_res.0.2)
lymphocyte$SCT_snn_res.0.3 <- as.character(lymphocyte$SCT_snn_res.0.3)
lymphocyte$SCT_snn_res.0.4 <- as.character(lymphocyte$SCT_snn_res.0.4)
lymphocyte$SCT_snn_res.0.5 <- as.character(lymphocyte$SCT_snn_res.0.5)
lymphocyte$SCT_snn_res.1 <- as.character(lymphocyte$SCT_snn_res.1)


## add prefixes
apcs$SCT_snn_res.0.1 <- paste0("A",apcs$SCT_snn_res.0.1)
apcs$SCT_snn_res.0.2 <- paste0("A",apcs$SCT_snn_res.0.2)
apcs$SCT_snn_res.0.3 <- paste0("A",apcs$SCT_snn_res.0.3)
apcs$SCT_snn_res.0.4 <- paste0("A",apcs$SCT_snn_res.0.4)
apcs$SCT_snn_res.0.5 <- paste0("A",apcs$SCT_snn_res.0.5)
apcs$SCT_snn_res.1 <- paste0("A",apcs$SCT_snn_res.1)

lymphocyte$SCT_snn_res.0.1 <- paste0("L",lymphocyte$SCT_snn_res.0.1)
lymphocyte$SCT_snn_res.0.2 <- paste0("L",lymphocyte$SCT_snn_res.0.2)
lymphocyte$SCT_snn_res.0.3 <- paste0("L",lymphocyte$SCT_snn_res.0.3)
lymphocyte$SCT_snn_res.0.4 <- paste0("L",lymphocyte$SCT_snn_res.0.4)
lymphocyte$SCT_snn_res.0.5 <- paste0("L",lymphocyte$SCT_snn_res.0.5)
lymphocyte$SCT_snn_res.1 <- paste0("L",lymphocyte$SCT_snn_res.1)

## get the mast cells barcodes
mast_annotation <- merged_samples$supercluster[merged_samples$supercluster == "Mast"] %>% names()


## save the objects
saveRDS(apcs, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/apcs_sub_sct.rds")
saveRDS(lymphocyte, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/lymphocyte_sub_sct.rds")


#<<<<<<<<<<<<<<<<<<<<<<< done LOCAL clustering



#>>>>>>>>>>>>>>>>>>>>>>> 
# global clustering
subcluster <- subset(merged_samples, subset = supercluster %in% c("APC","Lymphocyte","Mast")) %>% DietSeurat()

## run SCTransform and RunPCA on subclusters
subcluster <- SCTransform(subcluster, verbose=FALSE, method="glmGamPoi")
subcluster <- RunPCA(subcluster)

## run the scaledata, PCA for the combined subcluster.
subcluster <- subcluster %>% 
  ScaleData() %>% 
  RunPCA() 

subcluster %>% ElbowPlot(ndims = 50)

### Run harmony for the combined subcluster.
max_harmony_dims <- 20
subcluster <- subcluster %>% RunHarmony(group.by.vars = batch,
                                        dims.use = 1:max_harmony_dims, ncores=16)

ElbowPlot(subcluster, ndims = max_harmony_dims, reduction = 'harmony')

### recluster at resolutions 0.1, 0.2, 0.3, 0.4, 0.5 and 1
subcluster <- subcluster %>% recluster(dims = dims)


### add the annotation to the main "subcluster" object
global_0.1 <- subcluster$SCT_snn_res.0.1
global_0.2 <- subcluster$SCT_snn_res.0.2
global_0.3 <- subcluster$SCT_snn_res.0.3
global_0.4 <- subcluster$SCT_snn_res.0.4
global_0.5 <- subcluster$SCT_snn_res.0.5
global_1 <- subcluster$SCT_snn_res.1

mast_annotation <- merged_samples$supercluster[merged_samples$supercluster == "Mast"]


local_0.1_apcs <- apcs$SCT_snn_res.0.1
local_0.1_lymphocyte <- lymphocyte$SCT_snn_res.0.1
local_0.1 <- c(local_0.1_apcs, local_0.1_lymphocyte, mast_annotation)

local_0.2_apcs <- apcs$SCT_snn_res.0.2
local_0.2_lymphocyte <- lymphocyte$SCT_snn_res.0.2
local_0.2 <- c(local_0.2_apcs, local_0.2_lymphocyte, mast_annotation)

local_0.3_apcs <- apcs$SCT_snn_res.0.3
local_0.3_lymphocyte <- lymphocyte$SCT_snn_res.0.3
local_0.3 <- c(local_0.3_apcs, local_0.3_lymphocyte, mast_annotation)

local_0.4_apcs <- apcs$SCT_snn_res.0.4
local_0.4_lymphocyte <- lymphocyte$SCT_snn_res.0.4
local_0.4 <- c(local_0.4_apcs, local_0.4_lymphocyte, mast_annotation)

local_0.5_apcs <- apcs$SCT_snn_res.0.5
local_0.5_lymphocyte <- lymphocyte$SCT_snn_res.0.5
local_0.5 <- c(local_0.5_apcs, local_0.5_lymphocyte, mast_annotation)

local_1_apcs <- apcs$SCT_snn_res.1
local_1_lymphocyte <- lymphocyte$SCT_snn_res.1
local_1 <- c(local_1_apcs, local_1_lymphocyte, mast_annotation)


# add the metadata
subcluster <- AddMetaData(subcluster, metadata = global_0.5, col.name="Global_0.5")
subcluster <- AddMetaData(subcluster, metadata = global_0.4, col.name="Global_0.4")
subcluster <- AddMetaData(subcluster, metadata = global_0.3, col.name="Global_0.3")
subcluster <- AddMetaData(subcluster, metadata = global_0.2, col.name="Global_0.2")
subcluster <- AddMetaData(subcluster, metadata = global_0.1, col.name="Global_0.1")
subcluster <- AddMetaData(subcluster, metadata = global_1, col.name="Global_1")

subcluster <- AddMetaData(subcluster, metadata = local_0.5, col.name="Local_0.5")
subcluster <- AddMetaData(subcluster, metadata = local_0.4, col.name="Local_0.4")
subcluster <- AddMetaData(subcluster, metadata = local_0.3, col.name="Local_0.3")
subcluster <- AddMetaData(subcluster, metadata = local_0.2, col.name="Local_0.2")
subcluster <- AddMetaData(subcluster, metadata = local_0.1, col.name="Local_0.1")
subcluster <- AddMetaData(subcluster, metadata = local_1, col.name="Local_1")

# save the object
saveRDS(subcluster, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_sct.rds")



#>>>>>>>> Visualizations <<<<<<<<<<<#
# UMAPs
## General UMAP
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/UMAP_general.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="cluster", label=TRUE, raster=FALSE)+NoLegend()
dev.off()

## UMAPs for each custom resolution
clusters_labels <- c("Local_0.1","Local_0.2","Local_0.3","Local_0.4","Local_0.5","Local_1","Global_0.1","Global_0.2","Global_0.3","Global_0.4","Global_0.5","Global_1")

for(cl in clusters_labels){
  png(paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/UMAP_",cl,".png"), width=10, height=10, units="in", res=300)
  print(DimPlot(subcluster, group.by=cl, label=TRUE, raster=FALSE)+NoLegend())
  dev.off()
}


# Dot plots
## marker genes to generate dot plots
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


for(cl in clusters_labels){
  png(paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/dotplot_",cl,".png"), width=30, height=10, res=300, units="in")
  print(DotPlot(subcluster, features=rashx_markers, group.by=cl, cluster.idents=FALSE)+RotatedAxis())
  dev.off()
}


# Marker genes
## local
OUT <- createWorkbook()

for(cl in clusters_labels[1:6]){
  gc()
  Idents(subcluster) <- cl
  print(cl)
  markers <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "SCT", slot = "data", only.pos = TRUE)
  
  addWorksheet(OUT, paste0(cl))
	writeData(OUT, sheet = paste0(cl), x = markers)
}

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/marker_genes_local_resolutions.xlsx"))

## global
OUT <- createWorkbook()

for(cl in clusters_labels[7:12]){
  gc()
  Idents(subcluster) <- cl
  print(cl)
  markers <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "SCT", slot = "data", only.pos = TRUE)
  
  addWorksheet(OUT, paste0(cl))
	writeData(OUT, sheet = paste0(cl), x = markers)
}

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/marker_genes_global_resolutions.xlsx"))




# Count of cells per sample - local clusterign
OUT <- createWorkbook()
resolutions <- clusters_labels[1:6]

for(resolution in resolutions){
  print(resolution)
	Idents(subcluster) <- resolution
	my_list <- list()
	
	for(samplename in unique(subcluster$sample_id)){
		tmp <- subset(subcluster, subset = sample_id == samplename)
		my_list[[samplename]] <- as.data.frame(table(tmp[[resolution]]))
		colnames(my_list[[samplename]]) <- c("cluster","count_of_cells")
	}

## aggregate all the data frames in my_list into one data frame
	count_of_cells <- reduce(my_list, left_join, by="cluster")
	colnames(count_of_cells) <- c("cluster", names(my_list))

	addWorksheet(OUT, paste0(resolution))
	writeData(OUT, sheet = paste0(resolution), x = count_of_cells)
}

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/count_of_cells_per_sample_local.xlsx"))


# Count of cells per sample - global clusterign
OUT <- createWorkbook()
resolutions <- clusters_labels[7:12]

for(resolution in resolutions){
  print(resolution)
	Idents(subcluster) <- resolution
	my_list <- list()
	
	for(samplename in unique(subcluster$sample_id)){
		tmp <- subset(subcluster, subset = sample_id == samplename)
		my_list[[samplename]] <- as.data.frame(table(tmp[[resolution]]))
		colnames(my_list[[samplename]]) <- c("cluster","count_of_cells")
	}

## aggregate all the data frames in my_list into one data frame
	count_of_cells <- reduce(my_list, left_join, by="cluster")
	colnames(count_of_cells) <- c("cluster", names(my_list))

	addWorksheet(OUT, paste0(resolution))
	writeData(OUT, sheet = paste0(resolution), x = count_of_cells)
}

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global_sct/count_of_cells_per_sample_global.xlsx"))

# perform subclustering
merged_samples <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/merged_all_samples_integrated_and_annotated.RDS")
seed <- 123
set.seed(seed)

# variables to use in Harmony
batch <- c('condition', 'sample', 'chemistry', 'kit_10x')

# get some metadata
celltypes <- merged_samples@meta.data %>% dplyr::select(sample, sample_id, supercluster)

apcs_barcodes <- celltypes[celltypes$supercluster == "APC",] %>% rownames()
lymphocytes_barcodes <- celltypes[celltypes$supercluster == "Lymphocyte",] %>% rownames()

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# APCs and Lymphocytes LOCAL subclustering
apcs <- subset(merged_samples, subset = supercluster == "APC")
lymphocyte <- subset(merged_samples, subset = supercluster == "Lymphocyte")

## Findvariablegenes, scale and RunPCA
apcs <- apcs %>% DietSeurat() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()
lymphocyte <- lymphocyte %>% DietSeurat() %>% FindVariableFeatures(nfeatures = 3000) %>% ScaleData() %>% RunPCA()

## check the elbow plot
ElbowPlot(apcs)
ElbowPlot(lymphocyte)

## Run Harmony
max_harmony_dims <- 20
apcs <- apcs %>% RunHarmony(group.by.vars = batch,
                            dims.use = 1:max_harmony_dims, ncores=20)

apcs %>% ElbowPlot(ndims = max_harmony_dims, reduction = 'harmony')

lymphocyte <- lymphocyte %>% RunHarmony(group.by.vars = batch,
                            dims.use = 1:max_harmony_dims, ncores=20)

lymphocyte %>% ElbowPlot(ndims = max_harmony_dims, reduction = 'harmony')

## recluster after Harmony
dims <- 1:max_harmony_dims

### define the "recluster" function
recluster <- function(seuratobj, 
                      assay = 'RNA',
                      reduction = 'harmony',
                      dims = 1:20,
                      seed = 42,
                      n.neighbors = 30,
                      n.epochs = 200,
                      min.dist = 0.3,
                      k.param = n.neighbors,
                      metric = 'cosine',
                      algorithm = 1,
                      resolution = c(0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0),
                      default_ident = 'RNA_snn_res.1') {
  
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

## do the reclustering  at resolutions 0.5 and 1
apcs <- apcs %>% recluster(dims = dims, seed = seed, resolution = c(0.5,1))
lymphocyte <- lymphocyte %>% recluster(dims = dims, seed = seed, resolution = c(0.5,1))

## save the objects
saveRDS(apcs, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/apcs_sub.rds")
saveRDS(lymphocyte, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/lymphocyte_sub.rds")

## add prefixes
apcs$RNA_snn_res.0.5 <- paste0("A",apcs$RNA_snn_res.0.5)
apcs$RNA_snn_res.1 <- paste0("A",apcs$RNA_snn_res.1)

lymphocyte$RNA_snn_res.0.5 <- paste0("L",lymphocyte$RNA_snn_res.0.5)
lymphocyte$RNA_snn_res.1 <- paste0("L",lymphocyte$RNA_snn_res.1)

## add the sub-clustering @ 0.5 and 1 labels to the original object as new column
apcs_local_0.5_and_1 <- data.frame(local_0.5=apcs$RNA_snn_res.0.5, local_1=apcs$RNA_snn_res.1)
lymphocyte_local_0.5_and_1 <- data.frame(local_0.5=lymphocyte$RNA_snn_res.0.5, local_1=lymphocyte$RNA_snn_res.1)
B_annotation <- merged_samples$supercluster[merged_samples$supercluster == "B"] %>% names()
mast_annotation <- merged_samples$supercluster[merged_samples$supercluster == "Mast"] %>% names()

#<<<<<<<<<<<<<<<<<<<<<<< done LOCAL clustering



#>>>>>>>>>>>>>>>>>>>>>>> 
# global clustering
subcluster <- subset(merged_samples, subset = supercluster %in% c("APC","B","Lymphocyte","Mast")) %>% DietSeurat()

## assign the VariableFeatures for the subclusters.
combined_features <- union(apcs@assays$RNA@var.features, lymphocyte@assays$RNA@var.features)

VariableFeatures(subcluster) <- combined_features

## run the scaledata, PCA for the combined subcluster.
subcluster <- subcluster %>% 
  ScaleData() %>% 
  RunPCA() 

subcluster %>% ElbowPlot(ndims = 50)

### Run harmony for the combined subcluster.
max_harmony_dims <- 20
subcluster <- subcluster %>% RunHarmony(group.by.vars = batch,
                                        dims.use = 1:max_harmony_dims, ncores=20)

ElbowPlot(subcluster, ndims = max_harmony_dims, reduction = 'harmony')

### recluster at resolutions 0.5 and 1
subcluster <- subcluster %>% recluster(dims = dims, seed = seed, resolution = c(0.5,1), n.epochs = 500)

### add the annotation to the main "subcluster" object
global_0.5 <- subcluster$RNA_snn_res.0.5
global_1 <- subcluster$RNA_snn_res.1

B_annotation <- merged_samples$supercluster[merged_samples$supercluster == "B"]
mast_annotation <- merged_samples$supercluster[merged_samples$supercluster == "Mast"]

local_0.5_apcs <- apcs$RNA_snn_res.0.5
local_0.5_lymphocyte <- lymphocyte$RNA_snn_res.0.5
local_0.5 <- c(local_0.5_apcs, local_0.5_lymphocyte, B_annotation, mast_annotation)

local_1_apcs <- apcs$RNA_snn_res.1
local_1_lymphocyte <- lymphocyte$RNA_snn_res.1
local_1 <- c(local_1_apcs, local_1_lymphocyte, B_annotation, mast_annotation)

subcluster <- AddMetaData(subcluster, metadata = global_0.5, col.name="Global_0.5")
subcluster <- AddMetaData(subcluster, metadata = global_1, col.name="Global_1")
subcluster <- AddMetaData(subcluster, metadata = local_0.5, col.name="Local_0.5")
subcluster <- AddMetaData(subcluster, metadata = local_1, col.name="Local_1")

# save the "merged_samples" object excluding Fibroblast - including Local and Global clustering information
saveRDS(subcluster, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_with_local_and_global_0.5_and_1.rds")

#>>>>>>>> Visualizations <<<<<<<<<<<#
# UMAPs
## General UMAP
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/UMAP_general.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="supercluster", label=TRUE, raster=FALSE)+NoLegend()
dev.off()

## Local 0.5
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/UMAP_local_0.5.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="Local_0.5", label=TRUE, raster=FALSE)+NoLegend()
dev.off()

## Local 1
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/UMAP_local_1.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="Local_1", label=TRUE, raster=FALSE, repel=TRUE)+NoLegend()
dev.off()

## Global 0.5
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/UMAP_global_0.5.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="Global_0.5", label=TRUE, raster=FALSE)+NoLegend()
dev.off()

## Global 1
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/UMAP_global_1.png", width=10, height=10, units="in", res=300)
DimPlot(subcluster, group.by="Global_1", label=TRUE, raster=FALSE, repel=TRUE)+NoLegend()
dev.off()

# Dot plots
## marker genes to generate dot plots
rashx_markers <- c("CD3D",
                  "CCR7", "SELL", "KLF2","CD69", #Tn/Tcm/Tmm
                  "ITGAE","CXCR6", #Trm
                  "CD4","FOXP3", "TIGIT", "CTLA4", "IL2RA", #Treg
                  "CD8A","CD8B","GZMB", "NKG7", "CCL5", #CTL
                  "KLRB1", #NK/ILC
                  "GZMA", "GZMK", "KLRC1", "EOMES", "PRF1","KLRD1","GNLY", "CCR8", #NK
                  "IKZF3", "TBX21", #ILC1
                  "GATA3","MAF", "PTGDR2", "HPGDS", #ILC2
                  "RORC", "IL23R", "IL1R1", "KIT", "AHR", #ILC3
                  "TNFRSF4", "CD96","TNFRSF18","PRDM1","PDCD1","LAG3", #ac or ex
                  "BATF", "SNHG12", "ZFAS1", #Tet
                  "TRAT1","RORA","IL7R",  #not specific
                  "HLA-DRA","HLA-DRB1","CD83","IDO1", #APC
                  "CD207","EPCAM", #LC
                "CD68","CEBPB","FCER1G","C1QB","C1QC","CD163","NR4A1","NR4A2","KLF4","CD14","S100A9", #Mac
                  "MS4A7","LYZ","SERPINA1", "IL23A","IL1B", "CXCL3", #Mono
                  "CD1C", "CLEC10A","CLEC9A", "XCR1","FSCN1", "LAMP3",  #DC
                  "IGKC","JCHAIN","CD79A","MS4A1","CD19", "MS4A2", "IGHG1", "IGHG4", #B
                 "TPSB2", "TPSAB1", "CPA3",  "SIGLEC6", # Mast
                  "MKI67","TOP2A", "CENPF", "UBE2C", #cycling
                  "ITGA4","NCR1","IL17A","IL17F", #not specific
                    'PECAM1',   #stromal
                  'KRT1','KRT10','KRT14', #keratinocyte
                  'COL1A1','COL3A1','COL1A2',   #fibroblast
                  'TAGLN','ACTA2','RGS5',#smooth muscle cells
                  'KRT19','KRT7','CLDN10', #epithelial
                 'AQP1','CLDN5','PLVAP','EMCN'  #'SELE',   #Vasular endothelial cells
                  )

## make Dot Plot for the four different resolutions
png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/dotplot_local_0.5.png", width=30, height=10, res=300, units="in")
DotPlot(subcluster, features=rashx_markers, group.by="Local_0.5", cluster.idents=FALSE)+RotatedAxis()
dev.off()

png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/dotplot_local_1.png", width=30, height=10, res=300, units="in")
DotPlot(subcluster, features=rashx_markers, group.by="Local_1", cluster.idents=FALSE)+RotatedAxis()
dev.off()

png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/dotplot_global_0.5.png", width=30, height=10, res=300, units="in")
DotPlot(subcluster, features=rashx_markers, group.by="Global_0.5", cluster.idents=FALSE)+RotatedAxis()
dev.off()

png("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/dotplot_global_1.png", width=30, height=10, res=300, units="in")
DotPlot(subcluster, features=rashx_markers, group.by="Global_1", cluster.idents=FALSE)+RotatedAxis()
dev.off()


# Marker genes
Idents(subcluster) <- "Local_0.5"
markers_local_0.5 <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "RNA", slot = "data", only.pos = TRUE)
write.table(markers_local_0.5, "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/marker_genes_local_0.5.txt",
  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

Idents(subcluster) <- "Local_1"
markers_local_1 <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "RNA", slot = "data", only.pos = TRUE)
write.table(markers_local_1, "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/marker_genes_local_1.txt",
  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

Idents(subcluster) <- "Global_0.5"
markers_global_0.5 <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "RNA", slot = "data", only.pos = TRUE)
write.table(markers_global_0.5, "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/marker_genes_global_0.5.txt",
  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

Idents(subcluster) <- "Global_1"
markers_global_1 <- FindAllMarkers(subcluster, logfc.threshold = 0.25, verbose = TRUE, assay = "RNA", slot = "data", only.pos = TRUE)
write.table(markers_global_1, "/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/marker_genes_global_1.txt",
  col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)



# Count of cells per sample - local clusterign
OUT <- createWorkbook()
resolutions <- c("Local_0.5","Local_1")

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

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/count_of_cells_per_sample_local.xlsx"))


# Count of cells per sample - global clusterign
OUT <- createWorkbook()
resolutions <- c("Global_0.5","Global_1")

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

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/count_of_cells_per_sample_global.xlsx"))

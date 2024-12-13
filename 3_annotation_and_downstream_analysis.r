# clustering and downstream analyses - using local_1 and global_1 chosen by Jeff

## load in the data
subcluster <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_with_local_and_global_0.5_and_1.rds")

## add the custom annotation provided by Jeff
subcluster$Local_1_cluster <- subcluster$Local_1
subcluster$Local_1_supercluster <- subcluster$Local_1

subcluster$Global_1_cluster <- as.character(subcluster$Global_1)
subcluster$Global_1_supercluster <- as.character(subcluster$Global_1)

### local_1_cluster
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "B"] <- "B-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A15"] <- "B-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L2"] <- "CTL-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L11"] <- "CTL-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L12"] <- "CTL-3"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L17"] <- "CTL-4"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A4"] <- "DC-1 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A7"] <- "DC-2 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A3"] <- "DC-3 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A11"] <- "DC-4 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A12"] <- "DC-5 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A0"] <- "DC-6 moDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A13"] <- "DC-7"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A16"] <- "DC-8 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A18"] <- "DC-9 migDC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A8"] <- "LC"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A5"] <- "Mac-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A10"] <- "Mac-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A9"] <- "Mac-3"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A2"] <- "Mac-4"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A17"] <- "Mac-5 (LC)"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A14"] <- "Mac-6"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A19"] <- "Mac-7"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L15"] <- "Mac-8 (Stress)"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "Mast"] <- "Mast"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A1"] <- "Mono-1 (InfoMono)"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "A6"] <- "Mono-2 ?Mac?"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L16"] <- "NK"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L13"] <- "T(Activated/Stress)-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L19"] <- "T(Activated/Stress)-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L8"] <- "Tet"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L4"] <- "Tmm-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L5"] <- "Tmm-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L9"] <- "Tmm-3"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L3"] <- "Tn/Tcm"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L1"] <- "Treg-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L6"] <- "Treg-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L0"] <- "Trm-1"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L7"] <- "Trm-2"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L10"] <- "Trm-3"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L14"] <- "Trm-4"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L20"] <- "Trm-5 (NK?)"
subcluster$Local_1_cluster[subcluster$Local_1_cluster == "L18"] <- "Trm-6 (cycling)"

### local_1_supercluster
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "B"] <- "B"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A15"] <- "B"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L2"] <- "CTL"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L11"] <- "CTL"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L12"] <- "CTL"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L17"] <- "CTL"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A4"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A7"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A3"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A11"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A12"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A0"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A13"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A16"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A18"] <- "DC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A8"] <- "LC"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A5"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A10"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A9"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A2"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A17"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A14"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A19"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L15"] <- "Mac"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "Mast"] <- "Mast"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A1"] <- "Mono"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "A6"] <- "Mono"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L16"] <- "NK"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L13"] <- "T(Activated/Stress)"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L19"] <- "T(Activated/Stress)"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L8"] <- "Tet"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L4"] <- "Tmm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L5"] <- "Tmm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L9"] <- "Tmm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L3"] <- "Tn/Tcm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L1"] <- "Treg"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L6"] <- "Treg"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L0"] <- "Trm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L7"] <- "Trm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L10"] <- "Trm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L14"] <- "Trm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L20"] <- "Trm"
subcluster$Local_1_supercluster[subcluster$Local_1_supercluster == "L18"] <- "Trm"


### <------- global_1_cluster -------> ###
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "20"] <- "B"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "1"] <- "CTL-1"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "9"] <- "CTL-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "22"] <- "cycling (Trm/Treg)"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "7"] <- "DC-1 (migDC?)"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "15"] <- "DC-2 (migDC)"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "19"] <- "DC-3"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "18"] <- "LC"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "10"] <- "Mac-1"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "14"] <- "Mac-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "24"] <- "Mac-3?"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "17"] <- "Mac-4"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "26"] <- "Mac-5"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "8"] <- "Mast-1"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "23"] <- "Mast-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "21"] <- "Mono/Mac?"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "25"] <- "NK"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "16"] <- "T(activated/Stress)-1"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "11"] <- "T(activated/Stress)-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "12"] <- "Tet"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "4"] <- "Tmm?/Trm"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "3"] <- "Tn/Tcm"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "2"] <- "Treg-1"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "6"] <- "Treg-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "0"] <- "Trm-1 (Th17)"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "5"] <- "Trm-2"
subcluster$Global_1_cluster[subcluster$Global_1_cluster == "13"] <- "Trm-3?"

### local_1_supercluster
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "20"] <- "B"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "1"] <- "CTL"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "9"] <- "CTL"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "22"] <- "cycling (Trm/Treg)"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "7"] <- "DC"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "15"] <- "DC"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "19"] <- "DC"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "18"] <- "LC"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "10"] <- "Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "14"] <- "Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "24"] <- "Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "17"] <- "Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "26"] <- "Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "8"] <- "Mast"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "23"] <- "Mast"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "21"] <- "Mono/Mac"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "25"] <- "NK"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "16"] <- "T(activated/Stress)"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "11"] <- "T(activated/Stress)"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "12"] <- "Tet"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "4"] <- "Tmm/Trm"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "3"] <- "Tn/Tcm"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "2"] <- "Treg"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "6"] <- "Treg"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "0"] <- "Trm"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "5"] <- "Trm"
subcluster$Global_1_supercluster[subcluster$Global_1_supercluster == "13"] <- "Trm"

## save the object with annotation
saveRDS(subcluster, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_with_local_and_global_0.5_and_1_with_annotation.rds")


# Count of cells in each of the four systems
# Count of cells per sample - global clusterign
subcluster$Local_1_cluster <- factor(subcluster$Local_1_cluster)
subcluster$Local_1_supercluster <- factor(subcluster$Local_1_supercluster)
subcluster$Global_1_cluster <- factor(subcluster$Global_1_cluster)
subcluster$Global_1_supercluster <- factor(subcluster$Global_1_supercluster)

OUT <- createWorkbook()
resolutions <- c("Local_1_cluster","Local_1_supercluster","Global_1_cluster","Global_1_supercluster")

for(resolution in resolutions){
  print(resolution)
	Idents(subcluster) <- resolution
	my_list <- list()
	
	for(samplename in unique(subcluster$sample_id)){
		tmp <- subset(subcluster, subset = sample_id == samplename)
		my_list[[samplename]] <- as.data.frame(table(tmp[[resolution]]))
    my_list[[samplename]]$Freq <- paste0(my_list[[samplename]]$Freq," (",round(my_list[[samplename]]$Freq*100/sum(my_list[[samplename]]$Freq),2),"%)")
		colnames(my_list[[samplename]]) <- c("cluster","count_and_percentage_of_cells")
	}

## aggregate all the data frames in my_list into one data frame
	count_of_cells <- reduce(my_list, left_join, by="cluster")
	colnames(count_of_cells) <- c("cluster", names(my_list))

	addWorksheet(OUT, paste0(resolution))
	writeData(OUT, sheet = paste0(resolution), x = count_of_cells)
}

saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/count_of_cells_per_sample_local_global_cluster_and_supercluster.xlsx"))


# Stacked Bar plot to compare clusters' percentages
## add the required metadata
subcluster$supercondition <- subcluster$condition
subcluster$supercondition <- gsub("_Janssen\\d*\\s\\w*","",subcluster$supercondition) %>% gsub("tx","",.)
subcluster$supercondition %>% table()

subcluster$subcondition <- subcluster$condition
subcluster$subcondition <- gsub("_Janssen\\d*", "", subcluster$subcondition) %>% gsub("-tx","",.)

## fetch the data we need from the object
supercondition_local_1_cluster <- FetchData(subcluster, c("supercondition","Local_1_cluster"))
supercondition_local_1_supercluster <- FetchData(subcluster, c("supercondition","Local_1_supercluster"))
supercondition_global_1_cluster <- FetchData(subcluster, c("supercondition","Global_1_cluster"))
supercondition_global_1_supercluster <- FetchData(subcluster, c("supercondition","Global_1_supercluster"))

subcondition_local_1_cluster <- FetchData(subcluster, c("subcondition","Local_1_cluster"))
subcondition_local_1_supercluster <- FetchData(subcluster, c("subcondition","Local_1_supercluster"))
subcondition_global_1_cluster <- FetchData(subcluster, c("subcondition","Global_1_cluster"))
subcondition_global_1_supercluster <- FetchData(subcluster, c("subcondition","Global_1_supercluster"))



brks <- c(0, 0.25, 0.5, 0.75, 1)

### supercondition_local_1_cluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_supercondition_local_1_cluster.pdf", width=15, height=5)
ggplot(supercondition_local_1_cluster, aes(Local_1_cluster))+
  geom_bar(aes(fill = supercondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_supercondition_local_1_cluster.pdf", width=15, height=5)
ggplot(supercondition_local_1_cluster, aes(Local_1_cluster))+
  geom_bar(aes(fill = supercondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### supercondition_local_1_supercluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_supercondition_local_1_supercluster.pdf", width=15, height=5)
ggplot(supercondition_local_1_supercluster, aes(Local_1_supercluster))+
  geom_bar(aes(fill = supercondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_supercondition_local_1_supercluster.pdf", width=15, height=5)
ggplot(supercondition_local_1_supercluster, aes(Local_1_supercluster))+
  geom_bar(aes(fill = supercondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### supercondition_global_1_supercluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_supercondition_global_1_supercluster.pdf", width=15, height=5)
ggplot(supercondition_global_1_supercluster, aes(Global_1_supercluster))+
  geom_bar(aes(fill = supercondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_supercondition_global_1_supercluster.pdf", width=15, height=5)
ggplot(supercondition_global_1_supercluster, aes(Global_1_supercluster))+
  geom_bar(aes(fill = supercondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### supercondition_global_1_cluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_supercondition_global_1_cluster.pdf", width=15, height=5)
ggplot(supercondition_global_1_cluster, aes(Global_1_cluster))+
  geom_bar(aes(fill = supercondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_supercondition_global_1_cluster.pdf", width=15, height=5)
ggplot(supercondition_global_1_cluster, aes(Global_1_cluster))+
  geom_bar(aes(fill = supercondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()

# >>>>>>> END OF SUPER CONDITION <<<<<<<


# <<<<<<< START OF SUBCONDITION >>>>>>>>

### subcondition_local_1_cluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_subcondition_local_1_cluster.pdf", width=15, height=5)
ggplot(subcondition_local_1_cluster, aes(Local_1_cluster))+
  geom_bar(aes(fill = subcondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_subcondition_local_1_cluster.pdf", width=15, height=5)
ggplot(subcondition_local_1_cluster, aes(Local_1_cluster))+
  geom_bar(aes(fill = subcondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### subcondition_local_1_supercluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_subcondition_local_1_supercluster.pdf", width=15, height=5)
ggplot(subcondition_local_1_supercluster, aes(Local_1_supercluster))+
  geom_bar(aes(fill = subcondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_subcondition_local_1_supercluster.pdf", width=15, height=5)
ggplot(subcondition_local_1_supercluster, aes(Local_1_supercluster))+
  geom_bar(aes(fill = subcondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### subcondition_global_1_supercluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_subcondition_global_1_supercluster.pdf", width=15, height=5)
ggplot(subcondition_global_1_supercluster, aes(Global_1_supercluster))+
  geom_bar(aes(fill = subcondition), position = "fill")+theme_minimal()+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_subcondition_global_1_supercluster.pdf", width=15, height=5)
ggplot(subcondition_global_1_supercluster, aes(Global_1_supercluster))+
  geom_bar(aes(fill = subcondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()




### subcondition_global_1_cluster 
#### stacked bar plot PERCENTAGES
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_percentages_subcondition_global_1_cluster.pdf", width=15, height=5)
ggplot(subcondition_global_1_cluster, aes(Global_1_cluster))+
  geom_bar(aes(fill = subcondition), position = "fill")+
  scale_y_continuous(breaks = brks, labels = scales::percent(brks))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Percentage")
dev.off()

#### stacked bar plot COUNTS
pdf("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/stacked_bar_plot_count_subcondition_global_1_cluster.pdf", width=15, height=5)
ggplot(subcondition_global_1_cluster, aes(Global_1_cluster))+
  geom_bar(aes(fill = subcondition))+theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab("Count")
dev.off()

# >>>>>>>>>> END OF SUBCONDITION <<<<<<<<<<<<<<<<<<<#

# export the new object
saveRDS(subcluster, "/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_with_local_and_global_0.5_and_1_with_annotation_and_extra_metadata.rds")


#--------------------------------------------------------------------------------------------------------------------------------------------------------#

subcluster <- readRDS("/special_projects/55_scalp_psoriasis_analysis_jeff/RData/subcluster_with_local_and_global_0.5_and_1_with_annotation_and_extra_metadata.rds")

# DEGs
## PV vs Scalp Psoriasis - generate the DEGs genes for all resolutions
resolutions <- c("Local_1_cluster","Local_1_supercluster","Global_1_cluster","Global_1_supercluster")

for(res in resolutions){
  OUT <- createWorkbook()

	Idents(subcluster) <- res
  clusters <- levels(Idents(subcluster))
  print(res)
  
  for(cluster in clusters){
    print(cluster)
    tryCatch({
		  addWorksheet(OUT, cluster)
		  #Idents(merged_samples) <- res
		  
      markers <- FindMarkers(subcluster, only.pos=FALSE, logfc.threshold=0, min.pct=0.1, 
        group.by="supercondition", ident.1 = "Psoriasis Vulgaris", ident.2 = "Scalp psoriasis", subset.ident=cluster)
		  
      markers$gene <- rownames(markers)

      writeData(OUT, sheet = cluster, x = markers)
	  }, error=function(e){message(paste0("We have error!",res))})
  }
  saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/DEGs_",res,"_PV_vs_Scalp_psoriasis_per_cluster.xlsx"))
}
  

## Scalp Psoriasis Mid vs Scalp Psoriasis Pre - generate the DEGs genes for all resolutions
resolutions <- c("Local_1_cluster","Local_1_supercluster","Global_1_cluster","Global_1_supercluster")

for(res in resolutions){
  OUT <- createWorkbook()

	Idents(subcluster) <- res
  clusters <- levels(Idents(subcluster))
  print(res)
  
  for(cluster in clusters){
    print(cluster)
    tryCatch({
		  addWorksheet(OUT, cluster)
		  #Idents(merged_samples) <- res
		  
      markers <- FindMarkers(subcluster, only.pos=FALSE, logfc.threshold=0, min.pct=0.1, 
        group.by="subcondition", ident.1 = "Scalp psoriasis mid", ident.2 = "Scalp psoriasis pre", subset.ident=cluster)
		  
      markers$gene <- rownames(markers)

      writeData(OUT, sheet = cluster, x = markers)
	  }, error=function(e){message(paste0("We have error!",res))})
  }
  saveWorkbook(OUT, paste0("/special_projects/55_scalp_psoriasis_analysis_jeff/analysis_local_and_global/DEGs_",res,"_Scalp_psoriasis_Mid_vs_Scalp_psoriasis_Pre_per_cluster.xlsx"))
}


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

seurat_object_GSE18 <- read_rds("D:/Data_Sets/GSE18/harmonazied.Rds")

Idents(seurat_object_GSE18)<- seurat_object_GSE18$seurat_clusters
GSE18_Epithelial_N <- subset(seurat_object_GSE18, idents=c("5","18","3","13","6","16","20","19","14"))


# Code to Convert the Seurat Cluster to  Manual Clusters Name For Other Epithelail Cell 


cluster_numbers <- GSE18_Epithelial_N$seurat_clusters
names <-as.character(cluster_numbers)
for (i in 1:length(cluster_numbers)){
  if(cluster_numbers[i]== "5"|cluster_numbers[i]== "18"|cluster_numbers[i]== "3" |cluster_numbers[i]== "13")
  {
    names[i] <- "Epithelial Cell"
  }else if(cluster_numbers[i]== "14"  ){
    names[i] <- "Endothelial Cell"
  }else if(cluster_numbers[i]== "6" |cluster_numbers[i]== "20" |cluster_numbers[i]== "16"|cluster_numbers[i]== "19"){
    names[i] <- "Fibroblast"
  }
}

GSE18_Epithelial_N$Manual_Annotations <- names
GSE18_Epithelial_N@meta.data
GSE18_Epithelial_N_Manual_Annotated <-GSE18_Epithelial_N$Manual_Annotations
Idents(GSE18_Epithelial_N)<-GSE18_Epithelial_N$Manual_Annotations


#if you want to save the object use the below code
saveRDS(GSE18_Epithelial_N,"D:/Data_Sets/GSE18/New/Epithelial/GSE18_Epithelial_N_Manual_Annotated.Rds")


Idents(GSE18_Epithelial_N) <- GSE18_Epithelial_N$Manual_Annotations
DimPlot(GSE18_Epithelial_N, reduction = "umap", label = TRUE, repel = TRUE)


GSE18_Epithelial_Normalized  <- NormalizeData(GSE18_Epithelial_N)
GSE18_Epithelial_Normalized  <- FindVariableFeatures(GSE18_Epithelial_Normalized )
GSE18_Epithelial_Normalized  <- ScaleData(GSE18_Epithelial_Normalized )
GSE18_Epithelial_Normalized   <- RunPCA(GSE18_Epithelial_Normalized )
GSE18_Epithelial_Normalized   <- FindNeighbors(GSE18_Epithelial_Normalized  , dims = 1:30)
GSE18_Epithelial_Normalized   <- FindClusters(GSE18_Epithelial_Normalized  , resolution = 0.5)
GSE18_Epithelial_Normalized   <- RunUMAP(GSE18_Epithelial_Normalized  , dims = 1:30, n.neighbors = 50)

GSE18_Epithelial_Normalized@meta.data

DimPlot(GSE18_Epithelial_Normalized , reduction = "umap", group.by = "Manual_Annotations", split.by = 'condition')
DimPlot(GSE18_Epithelial_Normalized , reduction = "umap", group.by = "samples")

FeaturePlot(GSE18_Epithelial_Normalized, features = c ("EPCAM","PTPRC"), label = TRUE)  
FeaturePlot(GSE18_Epithelial_N, features = c ("EPCAM","PTPRC"), label = TRUE)

saveRDS(GSE18_Epithelial_Normalized ,"D:/Data_Sets/GSE18/New/Epithelial/GSE18_Epithelial_Normalized.Rds")


#####################################################################################

#############  Immune CELL SUBSET AND  ANNOTATION 

##############

#################################################################################




Idents(seurat_object_GSE18)<- seurat_object_GSE18$seurat_clusters
GSE18_Immnune_N <- subset(seurat_object_GSE18, idents=c("0","1","8","9","4","7","22","21","12","15"))


# Code to Convert the Seurat Cluster to  Manual Clusters Name For Other Epithelail Cell 


cluster_numbers <- GSE18_Immnune_N$seurat_clusters
names <-as.character(cluster_numbers)
for (i in 1:length(cluster_numbers)){
  if(cluster_numbers[i]== "0"|cluster_numbers[i]== "9")
  {
    names[i] <- "B Cell"
  }else if(cluster_numbers[i]== "1"|cluster_numbers[i]== "11"|cluster_numbers[i]== "8" ){
    names[i] <- "T Cell"
  }else if(cluster_numbers[i]== "4" ) {
    names[i] <- "NK Cell"
  }
  else if(cluster_numbers[i]== "7" |cluster_numbers[i]== "22" |cluster_numbers[i]== "21"){
    names[i] <- "Myeloid Cell"
  }else if(cluster_numbers[i]== "12" |cluster_numbers[i]== "15" ){
    names[i] <- "Dendritic Cell"
  }
  
}

GSE18_Immnune_N$Manual_Annotations <- names
GSE18_Immnune_N@meta.data
GSE18_Immnune_N_Manual_Annotated <-GSE18_Immnune_N$Manual_Annotations
Idents(GSE18_Immnune_N)<-GSE18_Immnune_N$Manual_Annotations


#if you want to save the object use the below code
saveRDS(GSE18_Immnune_N,"D:/Data_Sets/GSE18/New/Immune/GSE18_Immnune_N_Manual_Annotated.Rds")


Idents(GSE18_Immnune_N) <- GSE18_Immnune_N$Manual_Annotations
DimPlot(GSE18_Immnune_N, reduction = "umap", label = TRUE, repel = TRUE)


GSE18_Immune_Normalized   <- NormalizeData(GSE18_Immnune_N)
GSE18_Immune_Normalized   <- FindVariableFeatures(GSE18_Immune_Normalized  )
GSE18_Immune_Normalized   <- ScaleData(GSE18_Immune_Normalized  )
GSE18_Immune_Normalized    <- RunPCA(GSE18_Immune_Normalized  )
GSE18_Immune_Normalized    <- FindNeighbors(GSE18_Immune_Normalized   , dims = 1:30)
GSE18_Immune_Normalized    <- FindClusters(GSE18_Immune_Normalized   , resolution = 0.5)
GSE18_Immune_Normalized    <- RunUMAP(GSE18_Immune_Normalized   , dims = 1:30, n.neighbors = 50)



DimPlot(GSE18_Immune_Normalized  , reduction = "umap", group.by = "Manual_Annotations", split.by = 'new_condition')
DimPlot(GSE18_Immune_Normalized  , reduction = "umap", group.by = "samples")

FeaturePlot(seurat_object_GSE18, features = c ("EPCAM","PTPRC"), label = TRUE)  
FeaturePlot(GSE18_Immune_Normalized, features = c ("EPCAM","PTPRC"), label = TRUE)

saveRDS(GSE18_Immune_Normalized  ,"D:/Data_Sets/GSE18/New/Immune/GSE18_Immune_Normalized.Rds")


#################################################

############# IMMUNCE CELL SUB SET AND ANNOTATION 

################# After MANAUAL ANNOTATION 

################################################




Idents(GSE18_Immune_Normalized) <- GSE18_Immune_Normalized$Manual_Annotations
DimPlot(GSE18_Immune_Normalized ,reduction = "umap", label = TRUE, repel = TRUE)


Idents(GSE18_Immune_Normalized) <- GSE18_Immune_Normalized$Manual_Annotations
GSE11_B <- subset(GSE18_Immune_Normalized, idents=c("B Cell"))
GSE11_T <- subset(GSE18_Immune_Normalized, idents=c("T Cell"))
GSE11_NK <- subset(GSE18_Immune_Normalized, idents=c("NK Cell"))
GSE11_M <- subset(GSE18_Immune_Normalized, idents=c("Myeloid Cell"))
GSE11_D <- subset(GSE18_Immune_Normalized, idents=c("Dendritic Cell"))

GSE18_Immune_Normalized@meta.data$Manual_Annotations

DimPlot(GSE11_M ,reduction = "umap", label = TRUE, repel = TRUE)


saveRDS(GSE11_B,"D:/Data_Sets/GSE18/New/Immune/GSE11_B.Rds")
saveRDS(GSE11_T,"D:/Data_Sets/GSE18/New/Immune/GSE11_T.Rds")
saveRDS(GSE11_NK,"D:/Data_Sets/GSE18/New/Immune/GSE11_NK.Rds")
saveRDS(GSE11_M,"D:/Data_Sets/GSE18/New/Immune/GSE11_M.Rds")
saveRDS(GSE11_D,"D:/Data_Sets/GSE18/New/Immune/GSE11_D.Rds")






#################################################

#################  B CELL SUB SET AND ANNOTATION 

################# After MANAUAL ANNOTATION 

################################################



GSE11_B_Normalized   <- NormalizeData(GSE11_B)
GSE11_B_Normalized   <- FindVariableFeatures(GSE11_B_Normalized  )
GSE11_B_Normalized   <- ScaleData(GSE11_B_Normalized  )
GSE11_B_Normalized    <- RunPCA(GSE11_B_Normalized  )
GSE11_B_Normalized    <- FindNeighbors(GSE11_B_Normalized   , dims = 1:30)
GSE11_B_Normalized    <- FindClusters(GSE11_B_Normalized   , resolution = 0.5)
GSE11_B_Normalized    <- RunUMAP(GSE11_B_Normalized   , dims = 1:30, n.neighbors = 50)



DimPlot(GSE11_B_Normalized  , reduction = "umap", group.by = "Manual_Annotations", split.by = 'new_condition')
DimPlot(GSE11_B_Normalized  , reduction = "umap", group.by = "samples")



#################################################

#################  T CELL SUB SET AND ANNOTATION 

################# After MANAUAL ANNOTATION 

################################################


T_cell <- subset(GSE11_immune, idents=c("T Cell"))
DimPlot(GSE11_immune, reduction="umap", label=TRUE, repel=TRUE)

T_cell <- NormalizeData(T_cell)
T_cell <- FindVariableFeatures(T_cell)
T_cell <- ScaleData(T_cell )
T_cell <- RunPCA(T_cell )
T_cell <- FindNeighbors(T_cell , dims = 1:30)
T_cell  <- FindClusters(T_cell , resolution = 0.5)
T_cell  <- RunUMAP(T_cell , dims = 1:30, n.neighbors = 50)
DimPlot(T_cell, reduction="umap", label=TRUE, repel=TRUE)




#GSE11_T <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE11_T.Rds")


DimPlot(GSE11_T, reduction = "umap", label = TRUE, repel = TRUE)
Idents(GSE11_T) <- GSE11_T$Manual_Annotations

T_cELL_GSE18@meta.data
GSE11_T <-JoinLayers(GSE11_T)

GSE11_T_Normalized   <- NormalizeData(GSE11_T)
GSE11_T_Normalized   <- FindVariableFeatures(GSE11_T_Normalized)
GSE11_T_Normalized   <- ScaleData(GSE11_T_Normalized)
GSE11_T_Normalized    <- RunPCA(GSE11_T_Normalized)
GSE11_T_Normalized    <- FindNeighbors(GSE11_T_Normalized   , dims = 1:30)
GSE11_T_Normalized    <- FindClusters(GSE11_T_Normalized   , resolution = 0.5)
GSE11_T_Normalized    <- RunUMAP(GSE11_T_Normalized   , dims = 1:30, n.neighbors = 50)

DimPlot(GSE11_T_Normalized  , reduction = "umap", group.by = "Manual_Annotations", split.by = 'new_condition')
DimPlot(GSE11_T_Normalized  , reduction = "umap", group.by = "samples")
GSE11_immune <- JoinLayers(GSE11_immune) 



######################################################################################
























library(Seurat)


install.packages(c("arrow", "BiocManager", "bitops", "cachem", "cli", "colorspace", "corrplot", "curl", "data.table", "digest", "expm", "fastmap", "hdf5r", 
                   "htmltools", "httpuv", "lme4", "locfit", "matrixStats", "metap", "minqa", "mvtnorm", "nloptr", "openssl", "openxlsx", "parallelly", "polyclip", 
                   "promises", "Rcpp", "RCurl", "reticulate", "rlang", "rmarkdown", "robustbase", "RSpectra", "RSQLite", "Seurat", "SeuratObject", 
                   "spatstat.explore", "spatstat.geom", "spatstat.random", "spatstat.sparse", "spatstat.utils", "survival", "xfun", "yaml", "yulab.utils"))

install.packages(c("Seurat","SeureatData","SeuratObject","SeuratWrappers"))

base_pkgs <- c("package:stats", "package:graphics", "package:grDevices", 
               "package:utils", "package:datasets", "package:methods", 
               "package:base")


lapply(setdiff(search()[grepl("^package:", search())], base_pkgs), 
       detach, character.only = TRUE, unload = TRUE)

library(c("Seurat","SeureatData","SeuratObject","SeuratWrappers"))

library(Seurat)
library(SeuratData)
library(SeuratObject)
library(SeuratWrappers)

GSE11_immune <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE18_Immnune_N_Manual_Annotated.Rds")
library(tidyverse)
library(dplyr)


GSE11_immune <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE18_Immnune_N_Manual_Annotated.Rds")

GSE11_immune <- JoinLayers(GSE11_immune)
cell_calculations <- function(seurat_integrated,name,path){
  conditions <- unique(seurat_integrated$new_condition)
  final <- NULL
  for(i in 1:length(conditions)){
    print(i)
    print(conditions[i])
    Idents(seurat_integrated) <- seurat_integrated@meta.data$new_condition
    rorpos <- subset(seurat_integrated, idents=conditions[i])
    Idents(rorpos) <- rorpos@meta.data$Manual_Annotations
    n_cells_rorpos <- FetchData(rorpos, vars = "ident") %>%
      dplyr::count(ident) %>%
      tidyr::spread(ident, n)
    
    
    new <- as.data.frame( t(n_cells_rorpos))
    
    col<- rownames(new)
    col2 <- new[,1]
    
    new[,1] <- col
    new[,2] <- col2
    
    colnames(new) <- c("Cell_Type", "Counts")
    sumz <- sum(new$Counts)
    new$percentage <- (new$Counts/sumz)*100
    
    write.xlsx(new,paste0(path,name,"_",conditions[i],"_","counts",".xlsx"))
    #return(new)#
    pie(new$Counts, labels = new$Cell_Type, main = paste0("Pie Char of Cell Counts ",conditions[i]))
  }
  
}

GSE11_immune$new_condition

cell_calculations(GSE11_immune,name =  "GSE11_Immune",path = "D:/Data_Sets/GSE18/New/Immune/outputs/")

GSE11_immune$Manual_Annotations

library(openxlsx)

GSE11_epi <- read_rds("D:/Data_Sets/GSE18/New/Epithelial/GSE18_Epithelial_N_Manual_Annotated.Rds")

cell_calculations(GSE11_epi,name =  "GSE11_EPithelial",path = "D:/Data_Sets/GSE18/New/Epithelial/outputs/")

specific_labeled_degs_finder <- function(ROR2SeeuratObject, cluster,lfc=0.5,path,the_factor,merger=F){
  if(merger==T){
    ROR2SeeuratObject$joined <- paste0(the_factor,"_",ROR2SeeuratObject$condition)
    Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$joined
    cluster_cond <- unique(ROR2SeeuratObject$condition)
    cluster1 <- paste0(cluster_cond[1],"_",cluster_cond[2])
    print(unique(ROR2SeeuratObject$joined))
    print(cluster1)
    print(cluster2)
    print(cluster3)
    cluster2 <- paste0(cluster_cond[2],"_",cluster_cond[3])
    cluster3 <- paste0(cluster_cond[1],"_",cluster_cond[3])
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = lfc)
    marker1 <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster2, ident.2 = cluster3, min.pct = 0.25, logfc.threshold = lfc)
    marker2 <-FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster3, min.pct = 0.25, logfc.threshold = lfc)
  }else{
    ROR2SeeuratObject <- ROR2SeeuratObject
    Idents(ROR2SeeuratObject) <- the_factor
    cluster_cond <- unique(ROR2SeeuratObject$condition)
    cluster1 <- paste0(cluster_cond[1],"_",cluster_cond[2])
    cluster2 <- paste0(cluster_cond[2],"_",cluster_cond[3])
    cluster3 <- paste0(cluster_cond[1],"_",cluster_cond[3])
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = lfc,test.use = "MAST")
    marker1 <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster2, ident.2 = cluster3, min.pct = 0.25, logfc.threshold = lfc,test.use = "MAST")
    marker2 <-FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster3, min.pct = 0.25, logfc.threshold = lfc,test.use = "MAST")
  }
  filt <- function(markerz){
    markerz <- marker[markerz$p_val<0.05,]
    markerz$Names <- rownames(markerz)
    markerz <- marker[order(-markerz$avg_log2FC), ]
    return(markerz)
  }
  #marker <- marker[marker$p_val<0.05,]
  #marker$Names <- rownames(marker)
  #marker <- marker[order(-marker$avg_log2FC), ]
  marker <- filt(marker)
  marker1 <- filt(marker1)
  marker2 <- filt(marker2)
  write.csv(marker, paste0(path,as.character(cluster1),"_",as.character(cluster2),"-","conditional",".csv"),test.use = "MAST")
  write.csv(marker1, paste0(path,as.character(cluster2),"_",as.character(cluster3),"-","conditional",".csv"),test.use = "MAST")
  write.csv(marker2, paste0(path,as.character(cluster1),"_",as.character(cluster3),"-","conditional",".csv"),test.use = "MAST")
  return(marker)
}

specific_labeled_degs_finder
GSE11_immune$old_conditions <-GSE11_immune$condition
GSE11_immune$condition <- GSE11_immune$new_condition
specific_labeled_degs_finder(ROR2SeeuratObject = GSE11_immune,path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs",the_factor = GSE11_immune$Manual_Annotations, merger = T)
##############################################################################



specific_labeled_degs_finder <- function(ROR2SeeuratObject, cluster,lfc=0.5,path,the_factor,idn1,idn2,merger=F,conditions){
  if(merger==T){
    ROR2SeeuratObject$joined <- paste0(the_factor,"_",conditions)
    Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$joined
    cluster1 <- paste0(cluster,"_",idn1)
    cluster2 <- paste0(cluster,"_",idn2)
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1, ident.2 = cluster2, min.pct = 0.25, logfc.threshold = lfc)
  }else{
    ROR2SeeuratObject <- ROR2SeeuratObject
    Idents(ROR2SeeuratObject) <- the_factor
    cluster1 <- paste0(cluster,"_",idn1)
    cluster2 <- paste0(cluster,"_",idn2)
    marker <- FindMarkers(ROR2SeeuratObject, ident.1 = cluster1,ident.2 =cluster2, min.pct = 0.25, logfc.threshold = lfc)
  }
  marker <- marker[marker$p_val<0.05,]
  marker$Names <- rownames(marker)
  marker <- marker[order(-marker$avg_log2FC), ]
  write.csv(marker, paste0(path,as.character(cluster),"-",idn1,"-",idn2,".csv"))
  return(marker)
}


for (i in unique(GSE11_immune$Manual_Annotations)){
  specific_labeled_degs_finder(GSE11_immune, cluster = i, path = "D:/Data_Sets/GSE18/New/Immune/outputs/DEGs/",the_factor = GSE11_immune$Manual_Annotations, merger = T)
  
}
unique(paste0(GSE11_epi$Manual_Annotations,"_",GSE11_epi$condition))
GSE11_epi$condition <- GSE11_epi$new_condition
GSE11_epi <- JoinLayers(GSE11_epi)
for (i in unique(GSE11_epi$Manual_Annotations)){
  specific_labeled_degs_finder(GSE11_epi, cluster = i, path = "D:/Data_Sets/GSE18/New/Epithelial/outputs/DEGs/",the_factor = GSE11_epi$Manual_Annotations, merger = T)
  
}

#################################################

#################  NK CELL SUB SET AND ANNOTATION 

################# After MANAUAL ANNOTATION 

################################################

GSE11_NK <- read_rds("D:/Data_Sets/GSE18/New/Immune/GSE11_NK.Rds")


GSE11_NK_Normalized   <- NormalizeData(GSE11_NK)
GSE11_NK_Normalized   <- FindVariableFeatures(GSE11_NK_Normalized)
GSE11_NK_Normalized   <- ScaleData(GSE11_NK_Normalized)
GSE11_NK_Normalized    <- RunPCA(GSE11_NK_Normalized)
GSE11_NK_Normalized    <- FindNeighbors(GSE11_NK_Normalized   , dims = 1:30)
GSE11_NK_Normalized    <- FindClusters(GSE11_NK_Normalized   , resolution = 0.5)
GSE11_NK_Normalized    <- RunUMAP(GSE11_NK_Normalized   , dims = 1:30, n.neighbors = 50)



DimPlot(GSE11_T_Normalized  , reduction = "umap", group.by = "Manual_Annotations", split.by = 'new_condition')
DimPlot(GSE11_T_Normalized  , reduction = "umap", group.by = "samples")



#################################################

#################  Myeloid CELL SUB SET AND ANNOTATION  (M)

################# After MANAUAL ANNOTATION 

################################################


GSE11_M_Normalized   <- NormalizeData(GSE11_M)
GSE11_M_Normalized   <- FindVariableFeatures(GSE11_M_Normalized)
GSE11_M_Normalized   <- ScaleData(GSE11_M_Normalized)
GSE11_M_Normalized    <- RunPCA(GSE11_M_Normalized)
GSE11_M_Normalized    <- FindNeighbors(GSE11_M_Normalized   , dims = 1:30)
GSE11_M_Normalized    <- FindClusters(GSE11_M_Normalized   , resolution = 0.1)
GSE11_M_Normalized    <- RunUMAP(GSE11_M_Normalized   , dims = 1:30, n.neighbors = 50)



DimPlot(GSE11_M_Normalized  , reduction = "umap", group.by = "Manual_Annotations", split.by = 'new_condition')
DimPlot(GSE11_M_Normalized  , reduction = "umap", group.by = "seurat_clusters")






seurat_object_GSE18@meta.data

Idents(seurat_object_GSE18) <- seurat_object_GSE18$seurat_clusters
DimPlot(seurat_object_GSE18, reduction = "umap", label = TRUE, repel = TRUE)

DimPlot(seurat_object_GSE18, reduction = "umap", label = TRUE, repel = TRUE, split.by = "new_condition")

FeaturePlot(seurat_object_GSE18, features = c ("EPCAM","PTPRC"), label = TRUE)  
FeaturePlot(GSE18_Epithelial, features = c ("EPCAM","PTPRC"), label = TRUE)  

FeaturePlot(seurat_object_GSE18, features = c ("EPCAM","PTPRC","CAV1","KRT19","CRE"), label = TRUE)  

FeaturePlot(seurat_object_GSE18, features = c ("CD19", "CD1D", "CD21", "CD23", "CD5", "IGHD"), label = TRUE)  


FeaturePlot(NK_cells, features = c ("TXNIP", "IL7R", "RPL9P9", "ETS1"), label = TRUE)
FeaturePlot(seurat_object_GSE18, features = c ("TCL1A","CR2"), label = TRUE)
FeaturePlot(GSE11_B_SR, features = c ("CXCR4","CD69"), label = TRUE)


Idents(seurat_object_GSE18 ) <- seurat_object_GSE18 $seurat_clusters
DimPlot(seurat_object_GSE18 ,reduction = "umap", label = TRUE, repel = TRUE)

#####################################
Epithelail_Cells 

####################################
Idents(seurat_object_GSE18)<- seurat_object_GSE18$seurat_clusters
GSE18_Epithelial <- subset(seurat_object_GSE18, idents=c("5","18","3","13","6","16","20","19","14","12"))

GSE18_Epithelial <- NormalizeData(GSE18_Epithelial)
GSE18_Epithelial <- FindVariableFeatures(GSE18_Epithelial)
GSE18_Epithelial <- ScaleData(GSE18_Epithelial)
GSE18_Epithelial  <- RunPCA(GSE18_Epithelial)
GSE18_Epithelial <- FindNeighbors(GSE18_Epithelial, , dims = 1:35)
GSE18_Epithelial <- FindClusters(GSE18_Epithelial, resolution = 0.5)
GSE18_Epithelial  <- RunUMAP(GSE18_Epithelial , dims = 1:30, n.neighbors = 50)

DimPlot(GSE18_Epithelial, reduction = "umap", group.by = "seurat_clusters")
DimPlot(GSE18_Epithelial, reduction = "umap", group.by = "samples")

saveRDS(GSE18_Epithelial,"D:/Data_Sets/GSE18/GSE18_Epithelial.Rds")


#####################################
#Immune_Cells 

####################################

GSE18_Immune <- subset(seurat_object_GSE18, idents=c("4","1","11","8","10","2","17","0","7","9","22","15","21"))

GSE18_Immune <- NormalizeData(GSE18_Immune)
GSE18_Immune <- FindVariableFeatures(GSE18_Immune)
GSE18_Immune <- ScaleData(GSE18_Immune)
GSE18_Immune  <- RunPCA(GSE18_Immune)
GSE18_Immune  <- FindNeighbors(GSE18_Immune , dims = 1:30)
GSE18_Immune  <- FindClusters(GSE18_Immune , resolution = 0.5)
GSE18_Immune  <- RunUMAP(GSE18_Immune , dims = 1:30, n.neighbors = 50)



DimPlot(GSE18_Immune, reduction = "umap", group.by = "seurat_clusters")
DimPlot(GSE18_Immune, reduction = "umap", group.by = "samples")

saveRDS(GSE18_Immune,"D:/Data_Sets/GSE18/GSE18_Immune_harmonazied.Rds")



GSE18_Immune_Harmonized <- read_rds("D:/Data_Sets/GSE18/GSE18_Immune_harmonazied.Rds")




library(SingleCellExperiment)

GSE11_B <-JoinLayers(GSE11_B)

annotation <- SingleR(test = as.SingleCellExperiment(GSE11_B), ref = ref, labels = ref$label.main)
GSE11_B$Main_Annotations <- annotation$labels

ref <- celldex::DatabaseImmuneCellExpressionData()
results <- SingleR(test = as.SingleCellExperiment(GSE11_B ), ref = ref, labels = ref$label.main)
GSE11_B $singlr_labels <- results$labels
DimPlot(GSE11_B , reduction = 'umap', group.by = 'singlr_labels', label = TRUE, raster=FALSE)

Idents(GSE11_B) <-GSE11_B$seurat_clusters
Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$sub_annotations

GSE11_B@meta.data

Idents(GSE11_B) <- GSE11_B$singlr_labels
DimPlot(GSE11_B ,reduction = "umap", label = TRUE, split.by = 'condition')


GSE11_B <- subset(GSE11_B, idents=c("B cells"))
DimPlot(GSE11_B ,reduction = "umap", label = TRUE)


##################################################
#SingleR ANNOTATION

#################################################

GSE18_Epithelial_Harmonized_Corrected <- read_rds("D:/Data_Sets/GSE18/GSE18_Epithelial_harmonazied.Rds")

GSE18_Epithelial_Harmonized_Corrected <-JoinLayers(GSE18_Epithelial_Harmonized_Corrected)

annotation <- SingleR(test = as.SingleCellExperiment(GSE18_Epithelial_Harmonized_Corrected), ref = ref, labels = ref$label.main)
GSE18_Epithelial_Harmonized_Corrected$Main_Annotations <- annotation$labels

GSE18_Epithelial <-JoinLayers(GSE18_Epithelial)
ref <- celldex::HumanPrimaryCellAtlasData()
results <- SingleR(test = as.SingleCellExperiment(GSE18_Epithelial), ref = ref, labels = ref$label.main)
GSE18_Epithelial $singlr_labels <- results$labels
DimPlot(GSE18_Epithelial , reduction = 'umap', group.by = 'singlr_labels', label = TRUE, raster=FALSE)

M_CELL_GSE18 <-JoinLayers(M_CELL_GSE18)

ref <- celldex::DatabaseImmuneCellExpressionData()
results <- SingleR(test = as.SingleCellExperiment(GSE11_M_Monocytes_N ), ref = ref, labels = ref$label.main)
GSE11_M_Monocytes_N $singlr_labels <- results$labels
DimPlot(GSE11_M_Monocytes_N , reduction = 'umap', group.by = 'singlr_labels', label = TRUE, raster=FALSE)



#####################################################
##ANNOTATION
#####################################################

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
#ROR2SeeuratObject<-  read_rds("D:/Li-Porject Data/Li-Single Cell Data/outputs/Different_Resolutions/immuneResolution_9_30.Rds")

ROR2SeeuratObject <-NK_cells
#db_ = "D:/Li-Porject Data/Li-Single Cell Data/ROR2/Marker_List_Mouse_ROR2.xlsx"
#db_ = "D:/Li-Porject Data/Li-Single Cell Data/ROR2/Corrected/V2_final.xlsx"
#db_ = "D:/Li-Porject Data/Li-Single Cell Data/ROR2/Corrected/Updated Markers/V2_final_Paper.xlsx"
db_ = "D:/Sc_Type_Annotation/NK1.xlsx"
#names <- read.xlsx(db_)
#tissue = unique(names$tissueType)
tissue = c("Immune system")

gs_list = gene_sets_prepare(db_, tissue)


es.max = sctype_score(scRNAseqData = ROR2SeeuratObject[["RNA"]]$scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls = do.call("rbind", lapply(unique(ROR2SeeuratObject@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ROR2SeeuratObject@meta.data[ROR2SeeuratObject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ROR2SeeuratObject@meta.data$seurat_clusters==cl)), 10)
}))

# Top cell type for each cluster
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# Set low-confidence (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "NK_cell"
print(sctype_scores[,1:3])
#write.xlsx(sctype_scores,"D:/Li-Porject Data/Li-Single Cell Data/outputs/sctype_scores1.xlsx")

ROR2SeeuratObject@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  ROR2SeeuratObject@meta.data$sub_annotations[ROR2SeeuratObject@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(ROR2SeeuratObject, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sub_annotations', label.size = 5, raster = FALSE)
DimPlot(ROR2SeeuratObject, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)

DimPlot(ROR2SeeuratObject, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'seurat_clusters', label.size = 5)
#ggsave2("D:/Li-Porject Data/Li-Single Cell Data/outputs/ror/ROR_NEW.png")

ROR2SeeuratObject@meta.data

DimPlot(ROR2SeeuratObject, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sub_annotations', split.by= 'new_condition',label.size = 5)

########################
#if you want to save the object use the below code

#GSE11_B_Normalized <- ROR2SeeuratObject$sub_annotations
#Idents(GSE11_B_Normalized) <- ROR2SeeuratObject$sub_annotations
Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$sub_annotations
saveRDS(ROR2SeeuratObject,"D:/Data_Sets/GSE18/New/Immune/GSE11_NK_Normalized_ScType_Annotated.Rds")


#Myeloid data save 

Idents(ROR2SeeuratObject) <- ROR2SeeuratObject$sub_annotations
saveRDS(ROR2SeeuratObject,"D:/Data_Sets/GSE18/Immune/GSE11_M_Monocytes_ScType_Annotated.Rds")


############################## Check object created 
GSE11_B_Normalized_ScType_Annotated@meta.data

GSE11_B_Normalized_ScType_Annotated <- read_rds("D:/Data_Sets/GSE18/Immune/GSE11_B_Normalized_ScType_Annotated.Rds")

DimPlot(GSE11_B_Normalized_ScType_Annotated, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sub_annotations', split.by= 'condition',label.size = 5)
DimPlot(GSE11_B_Normalized_ScType_Annotated, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sub_annotations', split.by= 'new_condition',label.size = 5)


# DEGS OF CLUSTERS

NK_cells_Maker <- FindMarkers(NK_cells, ident.1 = 0, min.pct = 0.25 ) #lfc  ,logfc.threshold = 0.05
NK_cells_Maker <- NK_cells_Maker[NK_cells_Maker$p_val<0.05,]
NK_cells_Maker$Names <- rownames(NK_cells_Maker)
NK_cells_Maker <- NK_cells_Maker[order(-NK_cells_Maker$avg_log2FC), ]
write.csv(NK_cells_Maker, "D:/Data_Sets/GSE18/New/Immune/NK_Cluster0_New.csv")


#############################################
# Manual Annotation 

#########################################

##############################
#Immune_Cells 
#############################


# Code to Convert the Seurat Cluster to  Manual Clusters Name For Other Immune Cells  


cluster_numbers <- GSE18_Immune_Harmonized$seurat_clusters
names <-as.character(cluster_numbers)
for (i in 1:length(cluster_numbers)){
  if(cluster_numbers[i]== "11"|cluster_numbers[i]== "8"|cluster_numbers[i]== "14" |cluster_numbers[i]== "12")
  {
    names[i] <- "M"
  }else if(cluster_numbers[i]== "0" |cluster_numbers[i]== "3"|cluster_numbers[i]== "10" ){
    names[i] <- "B"
  }else if(cluster_numbers[i]== "4" |cluster_numbers[i]== "15" ){
    names[i] <- "NK"
  }else if(cluster_numbers[i]== "5" |cluster_numbers[i]== "7" |cluster_numbers[i]== "9"|cluster_numbers[i]== "13"|cluster_numbers[i]== "6"|cluster_numbers[i]== "16" |cluster_numbers[i]== "2"|cluster_numbers[i]== "1"){
    names[i] <- "T" 
  }
}

GSE18_Immune$Manual_Annotations <- names
GSE18_Immune@meta.data
GSE18_Immune_Manual_Annotated <-GSE18_Immune$Manual_Annotations
Idents(GSE18_Immune)<-GSE18_Immune$Manual_Annotations

#if you want to save the object use the below code
saveRDS(GSE18_Immune,"D:/Data_Sets/GSE18/Immune/GSE18_Immune_Manual_Annotated.Rds")

# Cluster No and Manaual_Annoataion Plot+ Condition Split 
DimPlot(GSE18_Immune,reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(GSE18_Immune,reduction = "umap", label = TRUE, group.by = 'seurat_clusters', repel = TRUE)
DimPlot(s.epi,reduction = "umap", label = TRUE, repel = TRUE,split.by = 'devs', group.by = 'seurat_clusters')

DimPlot(GSE18_Immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Manual_Annotations', label.size = 5)
DimPlot(GSE18_Immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'condition', label.size = 5)
DimPlot(GSE18_Immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Manual_Annotations',split.by = 'condition', label.size = 5)  
DimPlot(GSE18_Immune, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Manual_Annotations',split.by = 'samples', label.size = 5)  

GSE18_Immune@meta.data

##############################


# DEGS OF CLUSTERS

DimPlot(epi_3,reduction = "umap", label = TRUE, repel = TRUE)
epi_3Deg <- FindMarkers(epi_3, ident.1 = 1, min.pct = 0.25, logfc.threshold = 0.5)
epi_3Deg <- epi_3Deg[epi_3Deg$p_val<0.05,]
epi_3Deg$Names <- rownames(epi_3Deg)
epi_3Deg <- epi_3Deg[order(-epi_3Deg$avg_log2FC), ]
write.csv(epi_3Deg, "D:/Li-Porject Data/3-Mouse/Epi/cluster_1.csv")

# DEGS OF CLUSTERS Stroma 21
Com_All_DEG <- FindMarkers(Com_All, ident.1 = 21, min.pct = 0.25, logfc.threshold = 0.5)
Com_All_DEG <- Com_All_DEG[Com_All_DEG$p_val<0.05,]
Com_All_DEG$Names <- rownames(Com_All_DEG)
Com_All_DEG <- Com_All_DEG[order(-Com_All_DEG$avg_log2FC), ]
write.csv(Com_All_DEG, "D:/Li-Porject Data/3-Mouse/cluster_21_sTROMA_DEGS.csv")


############################
### GSE18 sUBSET&ANNOTATION
###########################D:/Data_Sets/GSE18/Immune/


GSE18_Immune_MA <- read_rds("D:/Data_Sets/GSE18/Immune/GSE18_Immune_Manual_Annotated.Rds")


Idents(GSE18_Immune_MA) <- GSE18_Immune_MA$Manual_Annotations
DimPlot(GSE18_Immune_MA ,reduction = "umap", label = TRUE, repel = TRUE)


#Idents(GSE19_B) <- GSE19_B$seurat_clusters
#DimPlot(GSE19_B, reduction = "umap", label = TRUE, repel = TRUE)
#GSE19_immune_updated <- subset(GSE19_immune, idents=c("0","1","2","3","4",
#"5","6","9"))

#Idents(GSE19_M) <-GSE19_M$seurat_clusters
Idents(GSE18_Immune_MA) <- GSE18_Immune_MA$Manual_Annotations
GSE11_B <- subset(GSE18_Immune_MA, idents=c("B"))
GSE11_T <- subset(GSE18_Immune_MA, idents=c("T"))
GSE11_NK <- subset(GSE18_Immune_MA, idents=c("NK"))
GSE11_M <- subset(GSE18_Immune_MA, idents=c("M"))
GSE11_T_NK <- subset(GSE18_Immune_MA, idents=c("T","NK"))

GSE18_Immune_MA@meta.data$Manual_Annotations

saveRDS(GSE11_B,"D:/Data_Sets/GSE18/Immune/GSE11_B.Rds")
saveRDS(GSE11_T,"D:/Data_Sets/GSE18/Immune/GSE11_T.Rds")
saveRDS(GSE11_NK,"D:/Data_Sets/GSE18/Immune/GSE11_NK.Rds")
saveRDS(GSE11_M,"D:/Data_Sets/GSE18/Immune/GSE11_M.Rds")


############################
### GSE18 sUBSET&ANNOTATION ### B cELL ANNOTATION AND Processing 
###########################D:/Data_Sets/GSE18/Immune/



GSE11_B@meta.data
DimPlot(GSE11_B, reduction = "umap", label = TRUE, repel = TRUE)
#T_basic_new_subtype <- subset(Immune_Basic,idents=c(1,2,3,4,9,11,14))

GSE11_B_SR <- NormalizeData(GSE11_B)
GSE11_B_SR <- FindVariableFeatures(GSE11_B_SR)
GSE11_B_SR <- ScaleData(GSE11_B_SR )
GSE11_B_SR  <- RunPCA(GSE11_B_SR )
GSE11_B_SR  <- FindNeighbors(GSE11_B_SR , dims = 1:30)
GSE11_B_SR  <- FindClusters(GSE11_B_SR , resolution = 0.3)
GSE11_B_SR  <- RunUMAP(GSE11_B_SR , dims = 1:30, n.neighbors = 50)




saveRDS(GSE11_B_SR,"D:/Data_Sets/GSE18/GSE11_B_Normalized_0.3.Rds")

GSE11_B_Normalized_0.3 <- read_rds("D:/Data_Sets/GSE18/GSE11_B_Normalized_0.3.Rds")
DimPlot(GSE11_B_Normalized_0.3, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'singlr_labels')
DimPlot(GSE11_B_Normalized_0.3, reduction = "umap", label = TRUE, repel = TRUE)
GSE11_B_Normalized@meta.data


############################
### GSE18 sUBSET&ANNOTATION ### T cELL ANNOTATION AND Processing 
###########################D:/Data_Sets/GSE18/Immune/

T_cELL_GSE18 <- read_rds("D:/Data_Sets/GSE18/Immune/GSE11_T.Rds")
DimPlot(T_cELL_GSE18, reduction = "umap", label = TRUE, repel = TRUE)
Idents(T_cELL_GSE18) <- T_cELL_GSE18$Manual_Annotations

T_cELL_GSE18@meta.data
GSE11_T <-JoinLayers(GSE11_T)
T_cELL_GSE18 <- NormalizeData(GSE11_T)
T_cELL_GSE18 <- FindVariableFeatures(T_cELL_GSE18)
T_cELL_GSE18 <- ScaleData(T_cELL_GSE18 )
T_cELL_GSE18  <- RunPCA(T_cELL_GSE18 )
T_cELL_GSE18  <- FindNeighbors(T_cELL_GSE18 , dims = 1:30)
T_cELL_GSE18  <- FindClusters(T_cELL_GSE18 , resolution = 0.5)
T_cELL_GSE18  <- RunUMAP(T_cELL_GSE18 , dims = 1:30, n.neighbors = 50)



############################
### GSE18 sUBSET&ANNOTATION ###  T_NK  CELL ANNOTATION AND Processing 
###########################D:/Data_Sets/GSE18/Immune/NK




GSE11_T <-JoinLayers(GSE11_T)
T_NK_CELL_GSE18 <- NormalizeData(GSE11_T_NK)
T_NK_CELL_GSE18 <- FindVariableFeatures(T_NK_CELL_GSE18)
T_NK_CELL_GSE18 <- ScaleData(T_NK_CELL_GSE18 )
T_NK_CELL_GSE18  <- RunPCA(T_NK_CELL_GSE18 )
T_NK_CELL_GSE18  <- FindNeighbors(T_NK_CELL_GSE18 , dims = 1:30)
T_NK_CELL_GSE18  <- FindClusters(T_NK_CELL_GSE18 , resolution = 0.5)
T_NK_CELL_GSE18  <- RunUMAP(T_NK_CELL_GSE18 , dims = 1:30, n.neighbors = 50)


############################
### GSE18 sUBSET&ANNOTATION ###  Myeloid  CELL ANNOTATION AND Processing 
###########################D:/Data_Sets/GSE18/Immune/Myeloid_Cells


M_CELL_GSE18 <- NormalizeData(GSE11_M)
M_CELL_GSE18 <- FindVariableFeatures(M_CELL_GSE18)
M_CELL_GSE18 <- ScaleData(M_CELL_GSE18 )
M_CELL_GSE18  <- RunPCA(M_CELL_GSE18 )
M_CELL_GSE18  <- FindNeighbors(M_CELL_GSE18 , dims = 1:30)
M_CELL_GSE18  <- FindClusters(M_CELL_GSE18 , resolution = 0.5)
M_CELL_GSE18  <- RunUMAP(M_CELL_GSE18 , dims = 1:30, n.neighbors = 50)


#MAKING OBJECT MORE CELAN BY SingleR Annotation and Cleaning


Idents(M_CELL_GSE18) <-M_CELL_GSE18$singlr_labels

GSE11_M_Monocytes <- subset(M_CELL_GSE18, idents=c("Monocytes"))
GSE11_M_Monocytes@meta.data

GSE11_M_Monocytes_N <- NormalizeData(GSE11_M_Monocytes)
GSE11_M_Monocytes_N <- FindVariableFeatures(GSE11_M_Monocytes_N)
GSE11_M_Monocytes_N <- ScaleData(GSE11_M_Monocytes_N )
GSE11_M_Monocytes_N  <- RunPCA(GSE11_M_Monocytes_N )
GSE11_M_Monocytes_N  <- FindNeighbors(GSE11_M_Monocytes_N , dims = 1:30)
GSE11_M_Monocytes_N  <- FindClusters(GSE11_M_Monocytes_N , resolution = 0.3)
GSE11_M_Monocytes_N  <- RunUMAP(GSE11_M_Monocytes_N , dims = 1:30, n.neighbors = 50)




M_CELL_GSE18@meta.data







samples <- seurat_object_GSE18$samples
new_samples <- samples
for (i in 1:length(samples)){
  if (grepl("cancer",samples[i])){
    new_samples[i] <- "Tumor"
  } else if(grepl("node 1",samples[i])){
    new_samples[i] <- "TLN 1"
  } else if(grepl("node 2",samples[i])){
    new_samples[i] <- "TLN 2"
  }
}

unique(new_samples)
seurat_object_GSE18$new_condition <- new_samples
saveRDS(seurat_object_GSE18,"D:/Data_Sets/GSE18/harmonazied.Rds")


Idents(seurat_object_GSE18)<- seurat_object_GSE18$seurat_clusters
GSE18_Epithelial <- subset(seurat_object_GSE18, idents=c("5","18","3","13","6","16","20","19","14","12"))
DimPlot(GSE18_Epithelial, reduction = "umap")
SE180286_Epithelial <- standard_processing_function(GSE18_Epithelial)


ref <- celldex::HumanPrimaryCellAtlasData()
SE180286_Epithelial <- JoinLayers(SE180286_Epithelial)
results <- SingleR(test = as.SingleCellExperiment(SE180286_Epithelial ), ref = ref, labels = ref$label.fine)

Idents(SE180286_Epithelial) <- SE180286_Epithelial$sun_annotations

DimPlot(SE180286_Epithelial, reduction = "umap")
GSE18_Immune <- subset(seurat_object_GSE18, idents=c("4","1","11","8","10","2","17","0","7","9","22","15","21"))

DimPlot(GSE18_Immune, reduction = "umap")
GSE18_Immune <- standard_processing_function(GSE18_Immune)


ref <- celldex::HumanPrimaryCellAtlasData()
GSE18_Immune <- JoinLayers(GSE18_Immune)
results <- SingleR(test = as.SingleCellExperiment(GSE18_Immune), ref = ref, labels = ref$label.fine)
GSE18_Immune$sub_annotations <- results$labels
Idents(GSE18_Immune) <- GSE18_Immune$sub_annotations

###############################################


install.packages("Signac")



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install



##################
GSE11_T_NK <- subset(GSE18_Immune_MA, idents=c("T","NK"))

GSE18_Immune_MA@meta.data$Manual_Annotations



GSE11_immune <- JoinLayers(GSE11_immune)



T_cell <- subset(GSE11_immune, idents=c("T Cell"))

T_cell <- NormalizeData(T_cell)
T_cell <- FindVariableFeatures(T_cell)
T_cell <- ScaleData(T_cell )
T_cell <- RunPCA(T_cell )
T_cell <- FindNeighbors(T_cell , dims = 1:30)
T_cell  <- FindClusters(T_cell , resolution = 0.3)
T_cell  <- RunUMAP(T_cell , dims = 1:30, n.neighbors = 50)
DimPlot(T_cell, reduction="umap", label=TRUE, repel=TRUE)


NK_cells <- subset(GSE11_immune, idents=c("NK Cell"))

NK_cells <- NormalizeData(NK_cells)
NK_cells <- FindVariableFeatures(NK_cells)
NK_cells <- ScaleData(NK_cells )
NK_cells <- RunPCA(NK_cells )
NK_cells <- FindNeighbors(NK_cells , dims = 1:30)
NK_cells  <- FindClusters(NK_cells , resolution = 0.3)
NK_cells  <- RunUMAP(NK_cells , dims = 1:30, n.neighbors = 50)
DimPlot(NK_cells, reduction="umap", label=TRUE, repel=TRUE)

T_cell
NK_cells

immune <- readRDS("D:/Data_Sets/GSE18/Immune/GSE18_Immune_Manual_Annotated.Rds")

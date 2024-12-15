###########################################################################################################
##THE AXUILARY FUNCTIONS
##########################################################################################################

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}










############################################################################################################################################################
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file
#
# @params: path_to_db_file - DB file with cell types
# @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#

gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}





#########################################################################################################################################################################
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# auto_detect_tissue_type: automatically detect a tissue type of the dataset
#
# @params: path_to_db_file - DB file with cell types
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells), 
# @params: scale - indicates whether the matrix is scaled (TRUE by default)

auto_detect_tissue_type <- function(path_to_db_file, scRNAseqData, scaled, ...){
  
  # get all tissue types in DB
  db_read = openxlsx::read.xlsx(path_to_db_file); tissues_ = unique(db_read$tissueType); result_ = c()
  
  for(tissue in tissues_){ print(paste0("Checking...", tissue));
    
    # prepare gene sets
    gs_list = gene_sets_prepare(path_to_db_file, tissue);
    
    es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = scaled, 
                          gs = gs_list$gs_positive, gs2 = gs_list$gs_negative, 
                          marker_sensitivity = gs_list$marker_sensitivity, verbose=!0);
    
    cL_resutls = do.call("rbind", lapply(unique(pbmc@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(pbmc@meta.data[pbmc@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    }))
    
    dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)
    
    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }
  
  # order by mean score
  result_ = result_[order(-result_$score),]
  
  # plot 
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",  main="The higher summary score, the more likely tissue type is")
  
  result_
}

##################################################################################################################
#END of Axuilary Functioins
##################################################################################################################

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
#ROR2SeeuratObject<-  read_rds("D:/Li-Porject Data/Li-Single Cell Data/outputs/Different_Resolutions/immuneResolution_9_30.Rds")

ROR2SeeuratObject <-s1
#db_ = "D:/Li-Porject Data/Li-Single Cell Data/ROR2/Marker_List_Mouse_ROR2.xlsx"
db_ = "D:/Li-Porject Data/Li-Single Cell Data/ROR2/Corrected/Updated Markers/V2_final_Paper_v1.xlsx"
tissue = c("Immune system","Epithelial Cell", "Endothelial Cell","Fibroblast")

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
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/5] = "Unknown"
print(sctype_scores[,1:3])
#write.xlsx(sctype_scores,"D:/Li-Porject Data/Li-Single Cell Data/outputs/sctype_scores1.xlsx")

ROR2SeeuratObject@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  ROR2SeeuratObject@meta.data$annotations[ROR2SeeuratObject@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(ROR2SeeuratObject, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotations', label.size = 5)
#ggsave2("D:/Li-Porject Data/Li-Single Cell Data/outputs/ror/ROR_NEW.png")


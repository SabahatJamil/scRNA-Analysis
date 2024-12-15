mini_kegg_OR <- function(instance_markers){
  instance_markers <- na.omit(instance_markers)
  rownames(instance_markers) <- instance_markers$Names
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  kegg_result <- enrichKEGG(gene = entrez_ids, organism = 'human', pvalueCutoff = 0.05)
  kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
  df <- as.data.frame(kegg_result@result)
  df <- df[df$p.adjust<0.05,]
  View(df)
}
mini_WP_OR <- function(instance_markers){
  instance_markers <- na.omit(instance_markers)
  rownames(instance_markers) <- instance_markers$Names
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  wp_result <- enrichWP(gene = entrez_ids, organism = "Homo sapiens", pvalueCutoff = 0.05)
  wp_result <- setReadable(wp_result, 'org.Hs.eg.db', 'ENTREZID')
  View(as.data.frame(wp_result@result))
  df <- as.data.frame(wp_result@result)
  df <- df[df$p.adjust<0.05,]
  View(df)
}
mini_Reac_OR <- function(instance_markers){
  instance_markers <- na.omit(instance_markers)
  rownames(instance_markers) <- instance_markers$Names
  instance_markers <- instance_markers[order(-instance_markers$avg_log2FC),]
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  reac_result <- enrichPathway(gene = entrez_ids, organism = "human", pvalueCutoff = 0.05)
  reac_result <- setReadable(reac_result, 'org.Hs.eg.db', 'ENTREZID')
  View(as.data.frame(reac_result@result))
  df <- as.data.frame(reac_result@result)
  df <- df[df$p.adjust<0.05,]
  View(df)
}

#enter the path to degs file csv
degs <- read_csv("D:/Data_Sets/Combined Object B and T/Annotated/Combined_B_T_DEGs/CD8_Cytotoxic-TLN-Tumor.csv")

#run the mini kegg for visualizing the results
mini_kegg_OR(degs)

#Enter the Description of each of the desired pathway
ctgry <- c("Th17 cell differentiation","Th1 and Th2 cell differentiation"
) 
#Cytokine-cytokine receptor interaction","Wnt signaling pathway","ECM-receptor interaction"

#this is the final step of this task to save the results
kegg_OR_GSEA(degs,name = "CD8_Cytotoxic-TLN-Tumor",path = "D:/Data_Sets/Combined Object B and T/Annotated/enrichment filtered/"
             ,ctgry =ctgry  )
unique(comb$sub_annotations)
Idents(comb) = comb$sub_annotations
B_cells <- subset(comb,idents=c("B_Plasma" ,     "B_Memory"   ,   "B_Activated"  , "B_Naive"  ,    "B_PB"))
b_degs$Names <- rownames(b_degs)
###similarly for wiki pathways
Idents(B_cells) = B_cells$updated_condtion
b_degs <- FindMarkers(B_cells,ident.1 = "TLN",ident.2 = "Tumor")
mini_kegg_OR(b_degs)
kegg_OR_GSEA(b_degs,name = "B_Cells_KEGG",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
             ,ctgry =c("T cell receptor signaling pathway","Th17 cell differentiation","Th1 and Th2 cell differentiation","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                       "FoxO signaling pathway","B cell receptor signaling pathway","TNF signaling pathway")  )

kegg_OR_GSEA(b_degs,name = "B_Cells_2",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
             ,ctgry =c("VEGF signaling pathway","p53 signaling pathway","NF-kappa B signaling pathway")  )

mini_WP_OR(b_degs)
WP_OR_GSEA(b_degs,name = "B_Cells",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("VEGFA-VEGFR2 signaling","	VEGFA-VEGFR2 signaling","T cell receptor and co-stimulatory signaling")  )
WP_OR_GSEA(b_degs,name = "B_Cells_2",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("IL2 signaling","IL4 signaling","IL7 signaling","IL1 signaling")  )
WP_OR_GSEA(b_degs,name = "B_Cells_3",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("Th17 cell differentiation pathway","IL11 signaling","IL18 signaling","IL6 signaling")  )




WP_OR_GSEA(b_degs,name = "B_Cells_Corrected",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("TNF-alpha signaling","B cell receptor signaling","T-cell receptor signaling","Interferon type I signaling",
                     "Modulators of TCR signaling and T cell activation","T cell receptor and co-stimulatory signaling","TNF-related weak inducer of apoptosis (TWEAK) signaling",
                     "IL2 signaling","IL4 signaling","IL7 signaling","IL1 signaling","IL11 signaling","IL18 signaling","IL6 signaling","Type II interferon signaling",
                     "IL3 signaling","IL7 signaling","	IL19 signaling")  )



mini_Reac_OR(b_degs)
Reactome_OR_GSEA(b_degs,name = "B_Cells_Reactome_OR",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
                 ,ctgry =c("Signaling by the B Cell Receptor (BCR)","Interferon Signaling","TNF signaling","	
Interleukin-17 signaling","Interleukin-1 signaling","Interferon alpha/beta signaling","PD-1 signaling","Interleukin-12 family signaling","Interleukin-3, Interleukin-5 and GM-CSF signaling") )





mini_WP_OR(degs)
WP_OR_GSEA(degs,name = "B_Cells",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("","") )


ctgry <- c("IL2 signaling","IL7 signaling","IL4 signaling","IL11 signaling","IL3 signaling","IL5 signaling","IL9 signaling","IL6 signaling",
           "Cancer immunotherapy by PD-1 blockade","Interferon type I signaling","TGF-beta receptor signaling","T-cell receptor signaling","B cell receptor signaling",
           "IL18 signaling","IL10 anti-inflammatory signaling","CCL18 signaling","TNF-alpha signaling","IL26 signaling","T cell receptor and co-stimulatory signaling",
           "IL19 signaling")


#imilarly for reactome
mini_Reac_OR(b_degs)
Reactome_OR_GSEA(b_degs,name = "B-Cell-Reactome",path = "D:/Data_Sets/Combined Object B and T/Annotated/enrichment filtered/"
                 ,ctgry =ctgry)

ctgry <- c("Interferon Signaling","TNF signaling","Downstream signaling events of B Cell Receptor (BCR)","Interleukin-17 signaling","Interleukin-1 signaling",
           "Signaling by Interleukins","Interferon alpha/beta signaling","	PD-1 signaling","Interleukin-12 signaling","Interleukin-3, Interleukin-5 and GM-CSF signaling")




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Idents(comb) = comb$sub_annotations
unique(comb$sub_annotations)
T_cells <- subset(comb,idents=c("CD4_Treg","CD4_Naive","CD4_Tfh","CD8_Activated","CD4_Tcm","CD8_Cytotoxic","CD8_Exhausted","CD8_Trm","CD4_Tem"))
Idents(T_cells) <- T_cells$updated_condtion
t_degs <- FindMarkers(T_cells,ident.1 = "TLN",ident.2 = "Tumor")


#GSEA is not going to be very very useful as it has only very limited amount of pathways to being with
#################################################################################################
instance_markers<- degs
gene_list <- instance_markers$avg_log2FC
names(gene_list) <- rownames(instance_markers)
gene_symbols <- instance_markers$Names
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
unamed_entrez_ids <- unname(entrez_ids)
newglist <- setNames(instance_markers$avg_log2FC, unamed_entrez_ids)
enrich_df <- gseKEGG(gene = newglist, organism = 'human', pvalueCutoff = 0.05)

enrich_df <- as.data.frame(enrich_df)
enrich_df <- enrich_df[enrich_df$p.adjust<0.05,]
#enrich_df <- enrich_df[enrich_df$Description %in% ctgry]
ggplot(enrich_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), 
                      size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_size_continuous(range = c(1, 4)) + 
  scale_color_gradientn(colors = c("blue", "yellow", "red")) + 
  theme_minimal() +
  labs(size = "Number", color = "-log10(Pvalue)", x = "GeneRatio", y = "Terms", 
       title = paste0(ident1," vs ",ident2)) +
  theme(axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"))

#################################################################################################

#T-Cell Processing for GSEA 


mini_kegg_OR(b_degs)
kegg_OR_GSEA(t_degs,name = "T_Cells_KEGG",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
             ,ctgry =c("T cell receptor signaling pathway","Th17 cell differentiation","Th1 and Th2 cell differentiation","PD-L1 expression and PD-1 checkpoint pathway in cancer",
                       "FoxO signaling pathway","B cell receptor signaling pathway","TNF signaling pathway")  )


mini_WP_OR(t_degs)
WP_OR_GSEA(b_degs,name = "B_Cells",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =c("VEGFA-VEGFR2 signaling","	VEGFA-VEGFR2 signaling","T cell receptor and co-stimulatory signaling")  )



library(DESeq2)
library(Seurat)




mini_kegg <- function(instance_markers){
  gene_list <- instance_markers$avg_log2FC
  names(gene_list) <- rownames(instance_markers)
  gene_symbols <- instance_markers$Names
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  unamed_entrez_ids <- unname(entrez_ids)
  newglist <- setNames(instance_markers$avg_log2FC, unamed_entrez_ids)
  # Perform KEGG pathway analysis
  kegg_result <- gseKEGG(gene = newglist, organism = 'human', pvalueCutoff = 0.05)
  
  View(as.data.frame(kegg_result@result) )
  
}
mini_kegg(degs)
#enter Description here in ctgry
ctgry <- c("","")
kegg_simp_plotter_human2(degs,name = "CD8_Cytotoxic-TLN-Tumor",path = "D:/Data_Sets/Combined Object B and T/Annotated/enrichment filtered/"
                         ,ctgry =ctgry)
###################################
#TESTING AREA
#Testing area
##############################################
condtion_degs <- FindMarkers(comb,ident.1 = "TLN", ident.2 = "Tumor")
condtion_degs$Names <- rownames(condtion_degs)
condtion_degs<- na.omit(condtion_degs)
condtion_degs <- condtion_degs[order(-condtion_degs$avg_log2FC),]
kegg_simp_plotter_human2(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
                         ,ctgry =10)
wp_simp_plotter_human2(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
                       ,ctgry =10)
rec_simp_plotter_human2(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
                        ,ctgry =10)
kegg_OR_GSEA(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
             ,ctgry =10)
WP_OR_GSEA(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
           ,ctgry =10)
Reactome_OR_GSEA(condtion_degs,name = "Overall",path = "D:/Data_Sets/Combined Object B and T/Annotated/Overall GSEA/"
                 ,ctgry =10)

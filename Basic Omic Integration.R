library(Signac)
library(Seurat)
library(GenomicRanges)
library(Matrix)

barcodes <- readLines(gzfile("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138_barcodes.tsv.gz"))
peaks <- read.table(gzfile("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138_peaks.bed.gz"), col.names = c("chr", "start", "end"))
peak_ranges <- makeGRangesFromDataFrame(peaks)
atac_data <- Read10X(data.dir = "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/")


matrix_file <- "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138_matrix.mtx.gz"
peaks_file <- "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138_peaks.bed.gz"
barcodes_file <- "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138_barcodes.tsv.gz"


peak_matrix <- Read10X(data.dir = "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/")


atac_data <- CreateSeuratObject(counts = peak_matrix, assay = "peaks")

####################################################
frag.file <- read.delim("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragments.tsv.gz", header = T)
counts <- Read10X("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/matrix.mtx")

data.dir <- "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/"
peak_matrix <- paste0(data.dir, "matrix.mtx.gz")
fragments <- paste0(data.dir, "fragments.tsv.gz")
barcodes <- paste0(data.dir, "barcodes.tsv.gz")
peaks1 <- paste0(data.dir, "peaks.bed.gz")

counts <- Read10X(data.dir = data.dir)
if (!requireNamespace("GenomicRanger", quietly = TRUE)) {
  BiocManager::install("GenomicRanger")
}
peaks <- read.table(peaks1, col.names = c("chr", "start", "end"))
gr_peaks <- makeGRangesFromDataFrame(peaks) 

peak_names <- paste(
  seqnames(gr_peaks), 
  start(gr_peaks), 
  end(gr_peaks), 
  sep = ":"
)

counts <- ReadMtx(
  mtx = peak_matrix,     
  features = peak_names,  
  cells = barcodes        
)

features_df <- data.frame(
  feature = peak_names,  
  stringsAsFactors = FALSE
)


features_file <- tempfile(fileext = ".txt")
write.table(features_df, file = features_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

counts <- ReadMtx(
  mtx = peak_matrix,
  features = features_file,
  cells = barcodes,
  feature.column = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,         
  sep = c(":", "-"),       
  genome = 'hg38',         
  fragments = fragments,   
  min.cells = 10           
)

counts <- ReadMtx(
  mtx = peak_matrix,
  features = gr_peaks,  
  cells = barcodes
)

if (!requireNamespace("ArchR", quietly = TRUE)) {
  BiocManager::install("ArchR")
}
library(ArchR)


arrow_path <- "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragment_index"


createArrowFile(
  inputFiles = fragments,
  outputDir = arrow_path,
  genome = "hg38",
  minFragments = 100
)
chrom_assay <- CreateChromatinAssay(
  counts = counts,         # Peak by cell count matrix
  sep = c(":", "-"),       # Separator for peak coordinates
  genome = 'hg38',         # Specify the genome (hg38 for human)
  fragments = fragments,   # Fragments file path
  min.cells = 10           # Minimum cells per peak
)


library(rtracklayer)

fragments <- read.table("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragments.tsv.gz", 
                        header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Create a GRanges object
gr_fragments <- GRanges(
  seqnames = fragments$V1,
  ranges = IRanges(start = fragments$V2, end = fragments$V3),
  score = fragments$V4
)

# Save the GRanges object as a RDS file
saveRDS(gr_fragments, "D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragment_index.rds")


fragment_index <- readRDS("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragment_index.rds")


chrom_assay <- CreateChromatinAssay(
  counts = counts,         
  sep = c(":", "-"),      
  genome = 'hg38',         
  fragments = fragment_index,  
  min.cells = 10           
)

fragment_index <- readRDS("D:/Data_Sets/GSE158398_atac/GSM4798906_P2003138/fragment_index.rds")


chrom_assay <- CreateChromatinAssay(
  counts = counts,         # Peak by cell count matrix
  sep = c(":", "-"),       # Separator for peak coordinates
  genome = 'hg38',         # Specify the genome (hg38 for human)
  fragments = fragment_index,  # Use the indexed fragment file
  min.cells = 10           # Minimum cells per peak
)


fragment_object <- CreateFragmentObject(
  path = fragments,  
  cells = barcodes                        
)



chrom_assay <- CreateChromatinAssay(
  counts = counts,        
  sep = c(":", "-"),      
  genome = 'hg38',        
  fragments = fragments,  
  min.cells = 10,         
  index.fragments = TRUE  
)





CreateChromatinAssay(index.f)




#######################################################################
folder_path <- "D:/Data_Sets/atac 2/"


file_names <- list.files(path = folder_path, full.names = TRUE)


print(file_names)

matrix_file <- "D:/Data_Sets/atac 2/GSE169246_TNBC_ATAC.counts.mtx.gz"
features_file <- "D:/Data_Sets/atac 2/GSE169246_TNBC_ATAC.feature.tsv.gz"
barcodes_file <- "D:/Data_Sets/atac 2/GSE169246_TNBC_ATAC.barcode.tsv.gz"

matrix_data <- ReadMtx(
  mtx = matrix_file,
  features = features_file,
  cells = barcodes_file
)


features_df <- read.table(features_file, header = FALSE, sep = "\t")
print(head(features_df))

matrix_data <- ReadMtx(
  mtx = matrix_file,
  features = features_file,
  cells = barcodes_file,
  feature.column = 1
)


atac_seurat <- CreateSeuratObject(counts = matrix_data, assay = "ATAC")


features_df <- read.table(features_file, header = FALSE, col.names = c("chr", "start", "end"))

features_df <- read.table(features_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)



chrom_assay <- CreateChromatinAssay(
  counts = matrix_data,
  sep = c(":", "-"),
  genome = 'hg38',         # Replace with your genome version if different
  min.cells = 10
)


atac_seurat[["ATAC"]] <- chrom_assay
View(atac_seurat@meta.data)

identical(rownames(chrom_assay), rownames(atac_seurat@assays$ATAC@counts))
identical(colnames(chrom_assay), colnames(atac_seurat@assays$ATAC@counts))



atac_seurat <- RunTFIDF(atac_seurat)
atac_seurat <- FindTopFeatures(atac_seurat, min.cutoff = 'q0')
atac_seurat <- RunSVD(atac_seurat)
atac_seurat <- RunUMAP(atac_seurat, reduction = 'lsi', dims = 2:30)
atac_seurat <- RunTSNE(atac_seurat, reduction = 'lsi', dims = 2:30)

atac_seurat <- FindNeighbors(atac_seurat, reduction = 'lsi', dims = 2:30)
atac_seurat <- FindClusters(atac_seurat, resolution = 0.5)
DimPlot(atac_seurat, reduction = "umap", label = TRUE)

saveRDS(atac_seurat ,"D:/Data_Sets/atac 2/atac_obj.Rds")
DimPlot(atac_seurat, reduction = "tsne", label = TRUE)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DirichletMultinomial")

library(chromVAR)
library(DirichletMultinomial)

BiocManager::install("JASPAR2020")
library(JASPAR2020)

atac_seurat <- readRDS("D:/Data_Sets/atac 2/atac_obj.Rds")

pwm_set <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# Add the motifs to the Seurat object
atac_seurat <- AddMotifs(
  object = atac_seurat,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pwm_set
)

atac_seurat <- RunChromVAR(atac_seurat, genome = BSgenome.Hsapiens.UCSC.hg38)


gene_activities <- GeneActivity(atac_seurat)
atac_seurat[['ACTIVITY']] <- CreateAssayObject(counts = gene_activities)
atac_seurat <- NormalizeData(atac_seurat, assay = 'ACTIVITY')
atac_seurat <- ScaleData(atac_seurat, assay = 'ACTIVITY')
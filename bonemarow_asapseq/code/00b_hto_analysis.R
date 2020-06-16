library(data.table)
library(Seurat)
library(Matrix)
library(BuenColors)

import_kite_counts <- function(){
  mtx <- fread(paste0("../data/hto/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../data/hto/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../data/hto/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}
hto <- import_kite_counts()
cells <- gsub("-1", "", fread("../data/barcodes/step1_barcodes_from_ArchR.tsv")[[1]])
cmat <- hto[,colnames(hto) %in% cells]
pbmc.hashtag <- CreateSeuratObject(counts = cmat, assay = "HTO")
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
pbmc.hashtag <- ScaleData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.995)
table(pbmc.hashtag$HTO_classification.global)
Idents(pbmc.hashtag) <- "HTO_maxID"
pR <- RidgePlot(pbmc.hashtag, assay = "HTO", 
                features = rownames(pbmc.hashtag[["HTO"]])[1:6], ncol = 3, cols = jdb_palette("corona"))
pR
#ggsave(pR, file = "HTO_ridges.pdf", width = 10, height = 5)  

# Visualize
Idents(pbmc.hashtag) <- "HTO_classification.global"
# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
pbmc.hashtag <- RunTSNE(pbmc.hashtag, distance.matrix = hto.dist.mtx, perplexity = 100)

DimPlot(pbmc.hashtag)

df <- data.frame(cells = paste0(colnames(pbmc.hashtag), "-1"), classification = (pbmc.hashtag$HTO_classification.global))
write.table(df[df$classification != "Doublet",1, drop = TRUE], file = "../data/barcodes/step2_singlets_from_HTO.tsv",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

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
doublet_scores_archr <-  fread("../data/barcodes/step1_barcodes_from_ArchR.tsv")[[2]]
vec <- doublet_scores_archr; names(vec) <- cells

# Subset cells
cmat <- hto[,colnames(hto) %in% cells]
bm.hashtag <- CreateSeuratObject(counts = cmat, assay = "HTO")
bm.hashtag <- NormalizeData(bm.hashtag, assay = "HTO", normalization.method = "CLR")
bm.hashtag <- ScaleData(bm.hashtag, assay = "HTO", normalization.method = "CLR")

bm.hashtag <- HTODemux(bm.hashtag, assay = "HTO", positive.quantile = 0.995)
val <- vec[rownames(bm.hashtag@meta.data)]
ddf <- data.frame(
  barcode = rownames(bm.hashtag@meta.data),
  classify_hto = bm.hashtag@meta.data$HTO_classification.global,
  val
)
ggplot(ddf, aes(x = classify_hto, y = val)) +
  geom_boxplot()

df <- data.frame(cells = paste0(colnames(bm.hashtag), "-1"), classification = (bm.hashtag$HTO_classification.global))
write.table(df[df$classification != "Doublet",1, drop = TRUE], file = "../data/barcodes/step2_singlets_from_HTO.tsv",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)



#---------------------------------------
Idents(bm.hashtag) <- "HTO_maxID"
pR <- RidgePlot(bm.hashtag, assay = "HTO", 
                features = rownames(bm.hashtag[["HTO"]])[1:6], ncol = 3, cols = jdb_palette("corona"))
pR
#ggsave(pR, file = "HTO_ridges.pdf", width = 10, height = 5)  

# Visualize
Idents(bm.hashtag) <- "HTO_classification.global"
# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = bm.hashtag, assay = "HTO"))))
# Calculate tSNE embeddings with a distance matrix
bm.hashtag <- RunTSNE(bm.hashtag, distance.matrix = hto.dist.mtx, perplexity = 100)

DimPlot(bm.hashtag)

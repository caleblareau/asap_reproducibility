library(data.table)
library(dplyr)
library(ArchR)
addArchRGenome("hg38")

import_kite_counts <- function(path, library, gz = FALSE){
  if(gz){
    xx <- ".gz"
  } else {
    xx <- ""
  }
  mtx <- fread(paste0(path, "featurecounts",library,".mtx",xx), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".barcodes.txt",xx), header = FALSE)[[1]], "-1")
  colnames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".genes.txt",xx), header = FALSE)[[1]])
  pct_human <- (matx[,2])/rowSums(matx)
  return(data.frame(barcode = rownames(matx), data.matrix(matx)))
}

addArchRGenome("hg38")
proj <- ArchRProject(
  ArrowFiles = c("../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_LLL.arrow", "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_OMNI.arrow"),
  outputDirectory = "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_combined_out",
  copyArrows = FALSE
)

dim(proj@cellColData)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", saveModel = FALSE)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")



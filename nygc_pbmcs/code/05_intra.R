library(data.table)
library(dplyr)
library(ArchR)
library(Seurat)

addArchRGenome("hg38")
source("00_import_kite.R")

proj <- ArchRProject(
  ArrowFiles = c("../../../asap_large_data_files/nygc_pbmc/output/archr/PBMCs_Intra.arrow"),
  outputDirectory = "../../../asap_large_data_files/nygc_pbmc/output/archr/intra_out",
  copyArrows = FALSE
)

set.seed(1)
dim(proj@cellColData)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", sampleCellsPre = 12000)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2
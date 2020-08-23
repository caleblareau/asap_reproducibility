library(data.table)
library(dplyr)
library(ArchR)
library(Seurat)

addArchRGenome("hg38")
source("00_import_kite.R")

proj <- ArchRProject(
  ArrowFiles = c("../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_LLL.arrow", "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_OMNI.arrow"),
  outputDirectory = "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_combined_out",
  copyArrows = FALSE
)

set.seed(1)
dim(proj@cellColData)
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", sampleCellsPre = 12000)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
saveRDS(proj, file = "../../../asap_large_data_files/nygc_pbmc/output/archr/ArchR_Proj_TSA-TSB_LLL-OMNI.rds")

# Append citeseq columns
LLL_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_LLL/featurecounts/", "", FALSE, "TSA_", "AvsB_LLL#")
OMNI_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_OMNI/featurecounts/", "", FALSE, "TSA_", "AvsB_OMNI#")
LLL_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_LLL/featurecounts/", "", FALSE, "TSB_", "AvsB_LLL#")
OMNI_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_OMNI/featurecounts/", "", FALSE, "TSB_", "AvsB_OMNI#")

mat_A <- rbind(LLL_kite_A, OMNI_kite_A)[rownames(proj@cellColData),]
mat_B <- rbind(LLL_kite_B, OMNI_kite_B)[rownames(proj@cellColData),]
all_adt <- cbind(mat_A, mat_B)

# Process the ADT in Seurat
pbmc <- CreateSeuratObject(counts = all_adt, assay = "ADT")
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay="ADT")
adt_mat_CLR <- data.matrix(pbmc@assays$ADT@data)
archr_df <- data.frame(
  proj@embeddings@listData$UMAP@listData$df,
  proj@cellColData,
  barcode = rownames(proj)
)

save(archr_df, all_adt, adt_mat_CLR, file = "../output/AvsB_LLL-Omni.rda")


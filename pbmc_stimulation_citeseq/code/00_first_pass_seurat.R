library(Seurat)

library(Seurat)
library(sctransform)
library(dplyr)
library(BuenColors)
library(Matrix)
library(data.table)
library(cowplot)
options(future.globals.maxSize = 4000 * 1024^2)

# Import
import_scRNAseq <- function(dir_base){
  
  name = dir_base
  data.dir <- paste0("../data/rnaseq/", dir_base)
  raw <- Read10X(data.dir = data.dir); colnames(raw) <- paste0(name, "-", colnames(raw))
  
  # Filter for singlet genes adn non-mitos
  raw <- raw[!grepl("^MT", rownames(raw)),]
  raw <- CreateSeuratObject(counts = raw,  min.cells = 3, min.features = 200); raw$source <- name; 
  raw <- NormalizeData(raw)
  raw <- FindVariableFeatures(raw)
  raw
  
}

# Import
pbmc_ctrl <- import_scRNAseq("ctrl")
pbmc_stim <- import_scRNAseq("stim")

features <- SelectIntegrationFeatures(object.list = list(pbmc_ctrl, pbmc_stim))
pbmc.list <- lapply(X = list(pbmc_ctrl, pbmc_stim), FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = pbmc.list, reduction = "rpca", 
                                  dims = 1:30)
pbmc.integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

# Do dimension reduction
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)

pbmc.integrated <- RunPCA(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:30, verbose = FALSE)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30, verbose = FALSE)
pbmc.integrated <- FindClusters(pbmc.integrated, verbose = FALSE, resolution = 0.5)

p1 <- DimPlot(pbmc.integrated, reduction = "umap", group.by = "source")
p2 <- DimPlot(pbmc.integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

FeaturePlot(pbmc.integrated, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                          "CCL2", "PPBP"), min.cutoff = "q9")

FeaturePlot(pbmc.integrated, features = c( "IFNG", "TNF"), split.by = "source", min.cutoff = "q9",
            cols = c("grey", "red"))

pbmc.integrated$cluster.stim <- paste(Idents(pbmc.integrated), pbmc.integrated$source, sep = "_")
Idents(pbmc.integrated) <- "cluster.stim"

clust0_response <- FindMarkers(pbmc.integrated, ident.1 = "7_stim", ident.2 = "7_ctrl", verbose = FALSE)
head(clust0_response, n = 15)

plots <- VlnPlot(pbmc.integrated, features = c("TNF", "IFNG", "CD69", "IL2"), split.by = "source", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(BuenColors)
library(harmony)

# Process kite counts
import_kite_counts <- function(library, bc_file, tech, bio){
  
  # Import the goodies
  bcs <- substr(fread(bc_file, header = FALSE)[[1]], 1, 16)
  mtx <- fread(paste0(library,"_featurecounts/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  
  # Label features
  rownames(matx) <- fread(paste0(library,"_featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(library,"_featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]])
  maty <- t(matx)[,rownames(matx) %in% bcs]
  maty <- maty[,Matrix::rowSums(maty) >= 500 & Matrix::rowSums(maty) < 500000]
  colnames(maty) <- paste0(library, "_", colnames(maty))
  
  # Now create a seurat object
  raw <- CreateSeuratObject(counts = maty,  min.cells = 3, min.features = 10, assay = "ADT"); raw$tech <- tech; raw$bio <- bio
  raw <- NormalizeData(raw)
  raw <- FindVariableFeatures(raw)

  return(raw)
}

# Import counts 
asap_control <- import_kite_counts("ASAP_ctrl",  "../../CONTROL_v12_hg38-mtMask_FC14k/outs/filtered_tf_bc_matrix/barcodes.tsv", "ASAP", "noStim")
asap_stim <- import_kite_counts("ASAP_stim",  "../../STIM_v12_hg38-mtMask_FC14k/outs/filtered_tf_bc_matrix/barcodes.tsv","ASAP", "Stim")
cite_control <- import_kite_counts("CITE_ctrl",  "../../../../asap_paper/asap_reproducibility/pbmc_stimulation_citeseq/data/rnaseq/ctrl/barcodes.tsv.gz","CITE", "noStim")
cite_stim <- import_kite_counts("CITE_stim",  "../../../../asap_paper/asap_reproducibility/pbmc_stimulation_citeseq/data/rnaseq/stim/barcodes.tsv.gz","CITE", "Stim")

features <- SelectIntegrationFeatures(object.list = list(asap_control, asap_stim, cite_control, cite_stim))
pbmc.list <- lapply(X = list(asap_control, asap_stim, cite_control, cite_stim), FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# Run CCA integration
anchors <- FindIntegrationAnchors(object.list = pbmc.list, reduction = "cca", 
                                  dims = 1:10)
pbmc.integrated <- IntegrateData(anchorset = anchors, dims = 1:10)

# Do dimension reduction
pbmc.integrated <- NormalizeData(pbmc.integrated, assay = "ADT", normalization.method = "CLR")
pbmc.integrated <- ScaleData(pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunPCA(pbmc.integrated, verbose = FALSE)

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:10)

pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:10, verbose = FALSE)
pbmc.integrated <- FindClusters(pbmc.integrated, verbose = FALSE, resolution = 0.5)
DimPlot(pbmc.integrated, group.by = "seurat_clusters", label = TRUE)

saveRDS(pbmc.integrated, file = "8JUNE2020_ADTseurat_PBMCstim.rds")

plot_df <- data.frame(pbmc.integrated@reductions$umap@cell.embeddings,
                      pbmc.integrated@meta.data)

ggplot(shuf(plot_df), aes(x = UMAP_1, y = UMAP_2, color = bio)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("grey", "firebrick"))

ggplot(shuf(plot_df), aes(x = UMAP_1, y = UMAP_2, color = tech)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))

pbmc.integrated$cluster.stim <- paste(Idents(pbmc.integrated), pbmc.integrated$bio, sep = "_")
Idents(pbmc.integrated) <- "cluster.stim"

FindMarkers(pbmc.integrated, ident.1 = "4_Stim", ident.2 = "4_noStim", verbose = FALSE) %>% head(10)

FeaturePlot(pbmc.integrated, features = c("CD3-1", "CD69", "CD56(NCAM)", "CD274", "CD25", "CD71", "CD4-1", "CD8",
                                          "CD19", "CD16", "CD14", "HLA-DR"),  min.cutoff = "q05", max.cutoff = "q95")


pbmc <- merge(merge(merge(asap_control, asap_stim,  add.cell.ids = c("a", "b")), cite_control, add.cell.ids= c("C1", "C2")), cite_stim,add.cell.ids = c("D1", "D2"))

# Do dimension reduction
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunHarmony(pbmc, c("tech","bio"), assay.use = "ADT")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:10)

pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE, reduction = "harmony")
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5)
DimPlot(pbmc)
plot_df2 <- data.frame(pbmc@reductions$umap@cell.embeddings,
                       pbmc@meta.data)

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = bio)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("grey", "firebrick"))

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = tech)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))

pbmc@meta.data$log2_nCount_ADT <- log2(pbmc@meta.data$nCount_ADT)
FeaturePlot(pbmc, "log2_nCount_ADT")

library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(viridis)


# Now process with module scores
pbmc_lll <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")
pbmc_lll@meta.data$stim <- ifelse(substr(colnames(pbmc_lll), 18, 18) == 1, "Control", "Stim")

DefaultAssay(pbmc_lll) <- "RNA"
pbmc_lll <- pbmc_lll  %>% RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction="pca", dims = 1:30) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) 

pbmc_lll <-   FindClusters(pbmc_lll, resolution = 0.1)
DimPlot(pbmc_lll, reduction = "umap", group.by = "stim", pt.size = .1, shuffle = TRUE)
DimPlot(pbmc_lll, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, label = TRUE)

DefaultAssay(pbmc_lll) <- "ADT"
FeaturePlot(pbmc_lll, features = c("CD3-1", "CD278", "CD69", "CD19", "CD56(NCAM)Recombinant"), split.by = "stim") &
  scale_color_viridis()

mod_scores<- fread("../output/LLL_module_scores.csv")
pbmc_lll@meta.data$protein_module_score <- mod_scores[["protein_MS1"]]
pbmc_lll@meta.data$rna_module_score <- mod_scores[["rna_MS1"]]
pbmc_lll@meta.data$chromatin_module_score <- mod_scores[["peaks.MS"]]


tcell_df <- pbmc_lll@meta.data %>% dplyr::filter(seurat_clusters %in% c(0, 1, 2, 3, 6))

ggplot(tcell_df, aes(x = chromatin_module_score, y = rna_module_score, color = protein_module_score)) + 
  geom_point() + scale_color_viridis()

cor(tcell_df[,c("chromatin_module_score", "rna_module_score", "protein_module_score")])


FeaturePlot(pbmc_lll, features = c("chromatin_module_score", "rna_module_score", "protein_module_score"), split.by = "stim") &
  scale_color_viridis()

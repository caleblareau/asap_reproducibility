library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(BuenColors)

pbmc_lll <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")

make_six_plot <- function(protein_name, gene_name){
  set.seed(1)
  DefaultAssay(pbmc_lll) = "ADT"
  protein_plot <- FeaturePlot(pbmc_lll, features = protein_name, split.by = "stim", ncol = 1, cells = sample(1:dim(pbmc_lll)[2]),
                              min.cutoff = "q02", max.cutoff = "q98", reduction = "wnn.3.umap")
  protein_control <- protein_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  protein_stim <- protein_plot[[2]] + theme_void() + theme(legend.position = "none")+ scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  DefaultAssay(pbmc_lll) = "GA"
  atac_plot <- FeaturePlot(pbmc_lll, features = gene_name, split.by = "stim", ncol = 1,cells = sample(1:dim(pbmc_lll)[2]),
                           min.cutoff = "q02", max.cutoff = "q98", reduction = "wnn.3.umap")
  atac_control <- atac_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red")) + ggtitle(NULL)
  atac_stim <- atac_plot[[2]] + theme_void() + theme(legend.position = "none")+scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  DefaultAssay(pbmc_lll) = "RNA"
  rna_plot <- FeaturePlot(pbmc_lll, features = gene_name, split.by = "stim", ncol = 1,cells = sample(1:dim(pbmc_lll)[2]),
                          min.cutoff = "q02", max.cutoff = "q98", reduction = "wnn.3.umap")
  rna_control <- rna_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  rna_stim <- rna_plot[[2]] + theme_void() + theme(legend.position = "none")+scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  cowplot::ggsave2(
    cowplot::plot_grid(atac_control, rna_control, protein_control, atac_stim, rna_stim, protein_stim, ncol = 3, scale = 0.9),
    width = 9, height = 6, file = paste0("../plots/dog/", protein_name, "_multiome_dog.png"))
  
}

make_six_plot("CD138-1(Syndecan-1)","SDC1")

make_six_plot("CD52","CD52")

make_six_plot("CD45-1","PTPRC")
make_six_plot("CD45-2","PTPRC")
make_six_plot("CD45RA","PTPRC")
make_six_plot("CD45RO","PTPRC")
make_six_plot("CD45RB","PTPRC")

make_six_plot("CD366", "TNFRSF18")
make_six_plot("CD357", "HAVCR2")


make_six_plot("CD3-1", "CD3E")
make_six_plot("CD3-2", "CD3E")


make_six_plot("CD8", "CD8A")
make_six_plot("CD4-1", "CD4")
make_six_plot("CD19", "CD19")
make_six_plot("CD20", "MS4A1")
make_six_plot("CD56(NCAM)", "NCAM1")
make_six_plot("CD184", "CXCR4")
make_six_plot("CD71", "TFRC")
make_six_plot("CD278", "ICOS")
make_six_plot("CD3-1", "CD3E")
make_six_plot("CD279", "PDCD1")
make_six_plot("CD69", "CD69")
make_six_plot("CD28", "CD28")
make_six_plot("CD25", "IL2RA")
make_six_plot("CD134","TNFRSF4")
make_six_plot("CD184", "CXCR4")
make_six_plot("CD71", "TFRC")
make_six_plot("CD142","F3")
make_six_plot("CD9","CD9")
make_six_plot("LOX-1","OLR1")
make_six_plot("CD223", "LAG3")
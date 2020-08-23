library(Signac)
library(Seurat)
library(BuenColors)

import_kite_counts <- function(library){
  mtx <- fread(paste0("../tag_data_LL/",library,"/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../tag_data_LL/",library,"/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../tag_data_LL/",library,"/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

ASAPB <- import_kite_counts("ADT3_ASAP_B")
colnames(ASAPB) <- paste0(colnames(ASAPB), "-2")

# Import cells
cells <- fread("../../../../../asap_large_data_files/broad_experiment_pbmcs/cellranger/filtered_peak_bc_matrix/barcodes.tsv", header = FALSE)[[1]]
cells_B <- cells[grepl("-2", cells)]

ASPAB_ss <- ASAPB[1:9, cells_B]

pbmc <- CreateSeuratObject(counts = ASPAB_ss, assay = "ADT")
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay="ADT")

# Process
cc <- fread("../../../../../asap_large_data_files/broad_experiment_pbmcs/cellranger/analysis/tsne/2_components/projection.csv")
boo <- substr(cc$Barcode, 18, 18) == "2"
dm <- data.matrix(cc[boo,c(2,3)])
colnames(dm) <- c("tsne_1", "tsne_2")

dr <- CreateDimReducObject(
  embeddings = dm, key = "tsne_", assay = "ADT"
)
DefaultAssay(pbmc) <- "ADT"
pbmc[['tsne']] <- dr

saveRDS(pbmc, file = "output/Broad_ASAPB_processed_ADT_object.rds")

p1 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD11c"),
            min.cutoff = "q20", max.cutoff = "q95", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p2 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD19"),
            min.cutoff = "q50", max.cutoff = "q90", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p3 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD3"),
            min.cutoff = "q10", max.cutoff = "q90", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p4 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD4"),
            min.cutoff = "q20", max.cutoff = "q80", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p5 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD8"),
            min.cutoff = "q20", max.cutoff = "q95", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p6 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD56"),
            min.cutoff = "q50", max.cutoff = "q95", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p7 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD16"),
            min.cutoff = "q20", max.cutoff = "q95", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

p8 <- FeaturePlot(pbmc, dims = c(1,2),
            features = c("CD14"),
            min.cutoff = "q20", max.cutoff = "q90", ncol = 3, reduction = "tsne") +
  theme_void() + ggtitle("") + theme(legend.position = "none")

cowplot::ggsave2(plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 1), 
                 file = "plots/ADT_8_viz.png", width = 15, height = 2.4, dpi = 300)

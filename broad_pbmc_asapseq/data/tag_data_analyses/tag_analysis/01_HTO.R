library(data.table)
library(Seurat)

import_kite_counts <- function(library){
  mtx <- fread(paste0("../tag_data_LL/",library,"/featurecounts.mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../tag_data_LL/",library,"/featurecounts.barcodes.txt"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../tag_data_LL/",library,"/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

ASAPC <- import_kite_counts("HTO1_ASAP_C")
cells <- gsub("-1", "", fread("../ASAP_C_v12_hg38-mtMask/outs/filtered_tf_bc_matrix/barcodes.tsv")[[1]])
cmat <- ASAPC[10:13,cells]

pbmc.hashtag <- CreateSeuratObject(counts = cmat, assay = "HTO")
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.9995)
table(pbmc.hashtag$HTO_classification.global)

Idents(pbmc.hashtag) <- "HTO_maxID"
pR <- RidgePlot(pbmc.hashtag, assay = "HTO", 
                features = rownames(pbmc.hashtag[["HTO"]])[1:4], ncol = 4, cols = jdb_palette("corona"))
#ggsave(pR, file = "plots/HTO_ridges.pdf", width = 10, height = 2.5)  

FeatureScatter(pbmc.hashtag, feature1 = "hto_HTO-1", feature2 = "hto_HTO-2")

# Visualize
Idents(pbmc.hashtag) <- "HTO_classification.global"
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = pbmc.hashtag, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
pbmc.hashtag<- RunTSNE(pbmc.hashtag, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(pbmc.hashtag)
pbmc.hashtag <- ScaleData(pbmc.hashtag)
m <- t(pbmc.hashtag@assays$HTO@counts)
idx_mat <- model.matrix(~0 + as.factor(max.col(m)))
pbmc.hashtag@meta.data$pct_max_counts <- rowSums(idx_mat * m)/rowSums(m)
pbmc.hashtag@meta.data$pct_max_counts_g0.5 <- pbmc.hashtag@meta.data$pct_max_counts > 0.5
FeaturePlot(pbmc.hashtag, features = c("pct_max_counts_g0.5"))

pBase <- DimPlot(pbmc.hashtag.subset, pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") +
  scale_color_manual(values = c("red2", "black")) + 
  theme_void() + theme(legend.position = "none") + ggtitle("")

g3gsave(pBase, file = "plots/doublets_tSNE.png", dpi = 500, width = 2, height = 2)


p1 <- FeaturePlot(pbmc.hashtag.subset, "HTO-1", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")

p2 <- FeaturePlot(pbmc.hashtag.subset, "HTO-2", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")

p3 <- FeaturePlot(pbmc.hashtag.subset, "HTO-3", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")

p4 <- FeaturePlot(pbmc.hashtag.subset, "HTO-4", max.cutoff = 'q99', pt.size = 0.01) +
  labs(x = "t-SNE1", y = "t-SNE2") + theme_void() + theme(legend.position = "none") + ggtitle("")

ggsave(cowplot::plot_grid(p1, p2, p3, p4, nrow = 2), file = "plots/HTO_tSNE.png", dpi = 500, width = 4, height = 4)

p1 <- HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000) +
  scale_fill_gradientn(colors = jdb_palette("solar_extra")) +
  theme(legend.position = "none")
ggsave(p1, file = "plots/visualize_HTO_heatmap.png", width = 3, height = 1.8)

# Now visualize on the cellranger map
df <- data.frame(pbmc.hashtag@meta.data)
df$Barcode <- paste0(rownames(df), "-3")
#write.table(df, file = "hash_assignments.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
tsne <- fread("../../../../../asap_large_data_files/broad_experiment_pbmcs/cellranger/analysis/tsne/2_components/projection.csv")
mdf_tsne <- data.frame(merge(df, tsne, by = "Barcode"))

p1 <- ggplot(mdf_tsne, aes(x = TSNE.1, y = TSNE.2, color = hash.ID)) +
  geom_point_rast(size = 0.1, raster.dpi = 500) + pretty_plot(fontsize = 8) + theme_void() + labs(color = "") +
  scale_color_manual(values = c("black", jdb_palette("corona")[1:4], "grey"))

cowplot::ggsave2(p1, file = "plots/tsne_color_hash.pdf", width = 2.4, height = 1.5)

                     
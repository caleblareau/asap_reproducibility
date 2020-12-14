library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)

gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

cells <- fread("ASAP_BROAD_3/outs/filtered_peak_bc_matrix/barcodes.tsv")[[1]]

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = "ASAP_BROAD_3/outs/fragments.tsv.gz",
  features = genebodyandpromoter.coords,
  cells = cells,
  chunk = 20
)
# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

saveRDS(gene.activities, file = "ASAP_BROAD_aggr_29April2020.rds")

pbmc <- CreateSeuratObject(
  counts = gene.activities,
  assay = 'RNA',
  project = 'RNA',
  min.cells = 1)

pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

cc <- fread("ASAP_BROAD_3/outs/analysis/tsne/2_components/projection.csv")
dm <- data.matrix(cc[,c(2,3)])
colnames(dm) <- c("tsne_1", "tsne_2")
rownames(dm) <- cc[[1]]
dm <- dm[colnames(pbmc),]
dr <- CreateDimReducObject(
  embeddings = dm, key = "tsne_", assay = "RNA"
)
DefaultAssay(pbmc) <- "RNA"
pbmc[['tsne']] <- dr

FeaturePlot(
  object = pbmc,
  features = c('MS4A1'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 1
)

cdf <- data.frame(
  Barcode = colnames(pbmc),
  library = substr(colnames(pbmc), 18, 18),
  pbmc@reductions$tsne@cell.embeddings,
  t(pbmc@assays$RNA@data[c('MS4A1', "TREM1", "LEF1"),])
)

p1t <- ggplot(cdf %>% filter(library == "1"), aes(x = tsne_1, y = tsne_2, color = TREM1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

p2t <- ggplot(cdf %>% filter(library == "2"), aes(x = tsne_1, y = tsne_2, color = TREM1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

p3t <- ggplot(cdf %>% filter(library == "3"), aes(x = tsne_1, y = tsne_2, color = TREM1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

p1m <- ggplot(cdf %>% filter(library == "1"), aes(x = tsne_1, y = tsne_2, color = MS4A1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

p2m <- ggplot(cdf %>% filter(library == "2"), aes(x = tsne_1, y = tsne_2, color = MS4A1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

p3m <- ggplot(cdf %>% filter(library == "3"), aes(x = tsne_1, y = tsne_2, color = MS4A1)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "blue"))

ggsave(cowplot::plot_grid(p1t, p1m, p2t,p2m,  p3t,  p3m, ncol = 2), 
       file = "plots/gene_score_viz.png", width = 3.66, height = 5, dpi = 500)


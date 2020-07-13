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
  mtx <- fread(paste0(library,"featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  
  # Label features
  rownames(matx) <- fread(paste0(library,"/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(library,"/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  maty <- t(matx)[,rownames(matx) %in% bcs]
  prop_control <- Matrix::colSums(maty[grepl("sotypeCtrl",rownames(maty)),])/Matrix::colSums(maty)
  maty <- maty[,Matrix::colSums(maty) >= 500 &prop_control < 0.01]
  colnames(maty) <- paste0(library, "_", colnames(maty))
  
  
  # Now create a seurat object
  raw <- CreateSeuratObject(counts = maty,  min.cells = 3, min.features = 10, assay = "ADT"); raw$tech <- tech; raw$bio <- bio
  raw <- NormalizeData(raw, assay = "ADT", normalization.method = "CLR")
  raw <- ScaleData(raw, assay = "ADT")
  
  return(raw)
}


# Import counts 
asap_control <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_ctrl_ADT/",  "../../pbmc_stimulation_asapseq/data/cellranger/control_barcodes.tsv", "ASAP", "noStim")
asap_stim <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_stim_ADT/",  "../../pbmc_stimulation_asapseq/data/cellranger/stim_barcodes.tsv","ASAP", "Stim")
cite_control <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/ctrl_featurecounts/",  "../../pbmc_stimulation_citeseq/data/rnaseq/ctrl/barcodes.tsv.gz","CITE", "noStim")
cite_stim <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/stim_featurecounts/",  "../../pbmc_stimulation_citeseq/data/rnaseq/stim/barcodes.tsv.gz","CITE", "Stim")

pbmc <- merge(merge(merge(asap_control, asap_stim,  add.cell.ids = c("AC", "AS")), cite_control, add.cell.ids= c("", "CC")), cite_stim,add.cell.ids = c("", "CS"))

# Do dimension reduction
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 100, selection.method = "vst")

# Manually curage variable features 
df <- pbmc@assays$ADT@meta.features; df$tag <- rownames(df)
tag_features <- df %>% arrange(desc(vst.variance)) %>% pull(tag) %>% head(150)
pbmc@assays$ADT@var.features <- tag_features

pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunHarmony(pbmc, c("tech", "bio"), assay.use = "ADT")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:15)

#pbmc <- FindNeighbors(pbmc, dims = 1:10, verbose = FALSE, reduction = "harmony")
#pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 0.5)

plot_df2 <- data.frame(pbmc@reductions$umap@cell.embeddings,
                       t(pbmc@assays$ADT@scale.data),
                       pbmc@meta.data)

library(viridis)
ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = CD3.2)) +
  facet_grid(tech~bio) +
  geom_point(size = 0.4) + pretty_plot() + labs(color = "CD3 CLR", x = "UMAP1", y = "UMAP2") +
  scale_color_viridis()

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = CD20)) +
  facet_grid(tech~bio) +
  geom_point(size = 0.4) + pretty_plot() + labs(color = "CD20 CLR", x = "UMAP1", y = "UMAP2") +
  scale_color_viridis()

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = CD11c)) +
  facet_grid(tech~bio) +
  geom_point(size = 0.4) + pretty_plot() + labs(color = "CD11c CLR", x = "UMAP1", y = "UMAP2") +
  scale_color_viridis()


ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = CD25)) +
  facet_grid(tech~bio) +
  geom_point(size = 0.4) + pretty_plot() + labs(color = "CD69 CLR", x = "UMAP1", y = "UMAP2") +
  scale_color_viridis()

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = bio)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "Condition", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("grey", "firebrick"))

ggplot(shuf(plot_df2), aes(x = UMAP_1, y = UMAP_2, color = tech)) +
  geom_point(size = 0.4) + pretty_plot() + L_border() + labs(color = "", x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))




pbmc@meta.data$log2_nCount_ADT <- log2(pbmc@meta.data$nCount_ADT)
FeaturePlot(pbmc, "log2_nCount_ADT")

library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(harmony)

# Import scRNA-seq
import_scRNAseq_cite_stim <- function(dir_base){
  
  data.dir <- paste0("../../pbmc_stimulation_citeseq/data/rnaseq/", dir_base)
  raw <- Read10X(data.dir = data.dir)
  colnames(raw) <- paste0(substr(colnames(raw), 1, 16), "-1")
  
  # import scrublet results
  singlets <- fread(paste0("../../pbmc_stimulation_citeseq/data/rnaseq/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(score < 0.2) %>% pull(barcode) # the original called threshold seemed too conservative; this is a better estimate for these libraries
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), singlets]
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- FindVariableFeatures(raw)
  raw <- NormalizeData(raw)
  raw <- ScaleData(raw)
  raw
}

# Import matrices
ctrl_scRNA <- import_scRNAseq_cite_stim("ctrl")
stim_scRNA <- import_scRNAseq_cite_stim("stim")
ga_control <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_control_ga.rds")
ga_stim <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_stim_ga.rds")

genes.use <- unique(c(VariableFeatures(ctrl_scRNA), VariableFeatures(stim_scRNA)))
length(genes.use)

seurat_integration_cca <- function(ga, pbmc.rna){
  
  # Process ATAC-seq data
  pbmc.atac <- CreateSeuratObject(counts = ga, assay = "ACTIVITY", project = "ATAC")
  pbmc.atac <- FindVariableFeatures(pbmc.atac)
  pbmc.atac <- NormalizeData(pbmc.atac)
  pbmc.atac <- ScaleData(pbmc.atac)
  
  # Establish transfer anchors
  pbmc.atac <- RunPCA(pbmc.atac, features = genes.use, verbose = FALSE)
  transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                          reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
  
  # Perform RNA imputation
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                             weight.reduction = pbmc.atac[["pca"]])
  pbmc.atac[["RNA"]] <- imputation
  coembed <- merge(x = pbmc.rna, y = pbmc.atac)
  
  # Now look at the co-embedding
  coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
  coembed <- RunUMAP(coembed, dims = 1:30)
  coembed
}
control_coembed <- seurat_integration_cca(ga_control, ctrl_scRNA)
stim_coembed <- seurat_integration_cca(ga_stim, stim_scRNA)

control_coembed$stim <- "control"
stim_coembed$stim <- "stim"

coembed4 <- merge(x = control_coembed, y = stim_coembed)
coembed4 <- ScaleData(coembed4, features = genes.use, do.scale = FALSE)
coembed4 <- RunPCA(coembed4, verbose = FALSE, features = genes.use)
coembed4 <- coembed4 %>%  RunHarmony(c("stim", "orig.ident"))
coembed4 <- coembed4 %>% RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(coembed4, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = c( 'stim'))

mdf <- data.frame(
  coembed4@reductions$umap@cell.embeddings,
  assay = coembed4@meta.data$orig.ident,
  stim = coembed4@meta.data$stim,
  cluster = coembed4@meta.data$seurat_clusters
)

set.seed(10)
ggplot(mdf,aes(x=UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point() + 
  pretty_plot() + 
  facet_grid(assay ~ stim) +
  scale_color_manual(values = sample(jdb_palette("corona")[1:16]))


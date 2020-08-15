library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(harmony)

source("10a_import_scRNA.R")

# Import matrices
ctrl_scRNA <- import_scRNAseq_cite_stim("ctrl"); dim(ctrl_scRNA)
stim_scRNA <- import_scRNAseq_cite_stim("stim"); dim(stim_scRNA)
ga_control <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_control_ga.rds")
ga_stim <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_stim_ga.rds")

# Process kite counts
import_kite_counts_known_bcs <- function(library, bcs, tech, bio){
  
  # Import the goodies
  bcs <- substr(bcs, 1, 16)
  mtx <- fread(paste0(library,"featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  
  # Label features
  rownames(matx) <- fread(paste0(library,"/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(library,"/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  maty <- t(matx)[,rownames(matx) %in% bcs]
  colnames(maty) <- paste0(colnames(maty), "-1")
  return(maty)
}

# Import counts 
cite_control <- import_kite_counts_known_bcs("../../pbmc_stimulation_citeseq/data/adt/ctrl_featurecounts/",  colnames(ctrl_scRNA),"CITE", "noStim")
cite_stim <- import_kite_counts_known_bcs("../../pbmc_stimulation_citeseq/data/adt/stim_featurecounts/",  colnames(stim_scRNA),"CITE", "Stim")
asap_control <- import_kite_counts_known_bcs("../../pbmc_stimulation_asapseq/data/adt/ASAP_ctrl_ADT/",  colnames(ga_control), "ASAP", "noStim")
asap_stim <- import_kite_counts_known_bcs("../../pbmc_stimulation_asapseq/data/adt/ASAP_stim_ADT/",  colnames(ga_stim),"ASAP", "Stim")

# Look for super expressors in the citeseq data
ctrl_meta <- data.frame(
  total = colSums(cite_control),
  control_tags = colSums(cite_control[grepl("sotypeCtrl", rownames(cite_control)),]) + 1
)

stim_meta <- data.frame(
  total = colSums(cite_stim),
  control_tags = colSums(cite_stim[grepl("sotypeCtrl", rownames(cite_stim)),]) + 1
)

ggplot(ctrl_meta, aes(total, control_tags)) +
  geom_point() + scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = 2.5e+04, color = "firebrick")+
  geom_hline(yintercept = 55, color = "firebrick")


ggplot(stim_meta, aes(total, control_tags)) +
  geom_point() + scale_x_log10() + scale_y_log10() +
  geom_vline(xintercept = 3e+04, color = "firebrick")+
  geom_hline(yintercept = 65, color = "firebrick")


remove_stim <- stim_meta$total > 25000 | stim_meta$control_tags > 55
remove_ctrl <- ctrl_meta$total > 30000 | ctrl_meta$control_tags > 65

# Now filter the citeseq data
ctrl_scRNA <- ctrl_scRNA[,!remove_ctrl]
stim_scRNA <- stim_scRNA[,!remove_stim]

cite_control <- cite_control[,!remove_ctrl]
cite_stim <- cite_stim[,!remove_stim]
dim(ctrl_scRNA); dim(cite_control)
dim(stim_scRNA); dim(cite_stim)

ctrl_scRNA[["ADT"]] <- CreateAssayObject(counts = cite_control)
stim_scRNA[["ADT"]] <- CreateAssayObject(counts = cite_stim)

# Figure out which genes to use in this analysis
gene_mapping <- fread("../data/marker_gene_mapping.tsv")
genes.use <- unique(c(VariableFeatures(ctrl_scRNA), VariableFeatures(stim_scRNA), gene_mapping[[2]]))
length(genes.use)
genes.use[genes.use %ni% rownames(ctrl_scRNA)]
genes.use <- genes.use[genes.use %in% rownames(ctrl_scRNA)]
length(genes.use)

seurat_integration_cca <- function(ga, pbmc.rna, protein){
  
  # Process ATAC-seq data
  pbmc.atac <- CreateSeuratObject(counts = ga, assay = "ACTIVITY", project = "ATAC")
  pbmc.atac[["ADT"]] <- CreateAssayObject(counts = protein)

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

control_coembed <- seurat_integration_cca(ga_control, ctrl_scRNA, asap_control)
stim_coembed <- seurat_integration_cca(ga_stim, stim_scRNA, asap_stim)

control_coembed$stim <- "control"
stim_coembed$stim <- "stim"

set.seed(1)
coembed4 <- merge(x = control_coembed, y = stim_coembed)
coembed4 <- ScaleData(coembed4, features = genes.use, do.scale = FALSE)
coembed4 <- RunPCA(coembed4, verbose = FALSE, features = genes.use)
coembed4@meta.data$ca <- paste0(coembed4@meta.data$orig.ident, "_", coembed4@meta.data$stim)
coembed4 <- NormalizeData(coembed4, assay = "ADT", normalization.method = "CLR")
coembed4 <- ScaleData(coembed4, assay = "ADT")
coembed4 <- coembed4 %>%  RunHarmony(c("stim", "orig.ident"))
coembed4 <- coembed4 %>%  FindNeighbors(reduction = "harmony", dims = 1:35)  %>% RunUMAP(reduction = "harmony", dims = 1:35)
coembed4 <- FindClusters(coembed4, resolution = 0.4) %>% identity()
DimPlot(coembed4, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = c( 'ca'), label = TRUE)

# Verify a couple of proteins
DefaultAssay(coembed4) = "ADT"
FeaturePlot(coembed4, reduction = "umap", features = c("CD4-1", "CD19", "CD16", "CD56(NCAM)", "CD8a", "CD14"), pt.size = .3, split.by = c('ca'))
saveRDS(coembed4, file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")


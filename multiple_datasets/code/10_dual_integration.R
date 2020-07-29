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

coembed4 <- merge(x = control_coembed, y = stim_coembed)
coembed4 <- ScaleData(coembed4, features = genes.use, do.scale = FALSE)
coembed4 <- RunPCA(coembed4, verbose = FALSE, features = genes.use)
coembed4 <- coembed4 %>%  RunHarmony(c("stim", "orig.ident"))
coembed4 <- coembed4 %>% RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(coembed4, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = c( 'stim'), label = TRUE)
FeaturePlot(coembed4, reduction = "umap", features = c("NCAM1", "CD3E", "MS4A1", "LYZ"), pt.size = .3, split.by = c( 'stim'))

mdf <- data.frame(
  coembed4@reductions$umap@cell.embeddings,
  assay = coembed4@meta.data$orig.ident,
  stim = coembed4@meta.data$stim,
  cluster = coembed4@meta.data$seurat_clusters,
  barcode = substr(rownames(coembed4@meta.data), 1, 16)
)

set.seed(10)
ggplot(mdf,aes(x=UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.1) + 
  pretty_plot() + 
  facet_grid(assay ~ stim) +
  scale_color_manual(values = sample(jdb_palette("corona")[1:16]))


t_cell_df <- mdf %>% filter(cluster %in% c(8,13)) # c(0, 1, 2, 3, 5, 6, 9, 10, 12)

# Quick function for counts per million
cpm <- function(mat){
  rs <- rowSums(mat)
  return(round(rs/sum(rs)*1000000, 2))
}

extract_cpm <- function(mat, barcodes){
  (cpm(mat[,substr(colnames(mat), 1, 16) %in% substr(barcodes, 1, 16)]))
}

# Assemble master table
atac_control <- extract_cpm(ga_control, t_cell_df %>% filter(assay == "ATAC" & stim == "control") %>% pull(barcode))
atac_stim <- extract_cpm(ga_stim, t_cell_df %>% filter(assay == "ATAC" & stim == "stim") %>% pull(barcode))

rna_control <- extract_cpm(ctrl_scRNA@assays$RNA@counts, t_cell_df %>% filter(assay == "RNA" & stim == "control") %>% pull(barcode))
rna_stim <- extract_cpm(stim_scRNA@assays$RNA@counts, t_cell_df %>% filter(assay == "RNA" & stim == "stim") %>% pull(barcode))

protein_control <- coembed4@assays$ADT@counts[, rownames(t_cell_df)[t_cell_df$stim == "control"]] %>% cpm
protein_stim <- coembed4@assays$ADT@counts[, rownames(t_cell_df)[t_cell_df$stim == "stim"]] %>% cpm

# Create changes for each assay
gene_mapping$protein_control <- protein_control[gene_mapping$Marker_name]
gene_mapping$protein_stim <- protein_stim[gene_mapping$Marker_name]
gene_mapping$log2_protein_change <- log2((gene_mapping$protein_stim + 1) / (gene_mapping$protein_control +  1))

gene_mapping$rna_control <- rna_control[gene_mapping$Gene_symbol]
gene_mapping$rna_stim <- rna_stim[gene_mapping$Gene_symbol]
gene_mapping$log2_rna_change <- log2((gene_mapping$rna_stim + 1) / (gene_mapping$rna_control +  1))

gene_mapping$atac_control <- atac_control[gene_mapping$Gene_symbol]
gene_mapping$atac_stim <- atac_stim[gene_mapping$Gene_symbol]
gene_mapping$log2_atac_change <- log2((gene_mapping$atac_stim + 1) / (gene_mapping$atac_control +  1))

p1 <- ggplot(gene_mapping %>% filter(rna_control > 10 | rna_stim > 10),
       aes(x = log2_rna_change, y = log2_protein_change, color = log2_atac_change, label = Marker_name)) +
  geom_point() + scale_color_gradientn(colors = jdb_palette("solar_basic")) +
  pretty_plot() + L_border() +theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(-5.5, 5.5)) + scale_y_continuous(limits = c(-5.5, 5.5))

p2 <- ggplot(gene_mapping %>% filter(rna_control > 10 | rna_stim > 10),
       aes(x = log2_rna_change, y = log2_protein_change,  label = Marker_name, color = log2_atac_change)) +
  geom_text(size = 3) + scale_color_gradientn(colors = jdb_palette("solar_basic")) +
  pretty_plot() + L_border() +theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(-5.5, 5.5)) + scale_y_continuous(limits = c(-5.5, 5.5))

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), file = "../plots/scatter_DOG_Myeloid.png", width = 8, height = 5)
write.table(gene_mapping, file = "../output/DOG_source_Myeloid.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




coembed4 <- NormalizeData(coembed4, assay = "ADT", normalization.method = "CLR")
coembed4 <- ScaleData(coembed4, assay = "ADT")

make_six_plot <- function(protein_name, gene_name){
  DefaultAssay(coembed4) = "ADT"
  protein_plot <- FeaturePlot(coembed4, features = protein_name, split.by = "stim", ncol = 1)
  protein_control <- protein_plot[[1]] + theme_void() + theme(legend.position = "none") + ggtitle("ADT control")
  protein_stim <- protein_plot[[2]] + theme_void() + theme(legend.position = "none")+ ggtitle("ADT stim")
  
  DefaultAssay(coembed4) = "ACTIVITY"
  coembedplot_atac_obj <- coembed4[,coembed4$orig.ident == "ATAC"]
  atac_plot <- FeaturePlot(coembedplot_atac_obj, features = gene_name, split.by = "stim", ncol = 1)
  atac_control <- atac_plot[[1]] + theme_void() + theme(legend.position = "none") + ggtitle("ATAC control")
  atac_stim <- atac_plot[[2]] + theme_void() + theme(legend.position = "none")+ ggtitle("ATAC stim")
  
  DefaultAssay(coembed4) = "RNA"
  coembedplot_rna_obj <- coembed4[,coembed4$orig.ident == "RNA"]
  rna_plot <- FeaturePlot(coembedplot_rna_obj, features = gene_name, split.by = "stim", ncol = 1)
  rna_control <- rna_plot[[1]] + theme_void() + theme(legend.position = "none") + ggtitle("RNA control")
  rna_stim <- rna_plot[[2]] + theme_void() + theme(legend.position = "none")+ ggtitle("RNA stim")
  
  cowplot::ggsave2(
    cowplot::plot_grid(atac_control, rna_control, protein_control, atac_stim, rna_stim, protein_stim, ncol = 3),
    width = 9, height = 6, file = paste0("../plots/dog/", protein_name, "_dog.png"))
  
}

all_dog <- gene_mapping[gene_mapping$Gene_symbol %in% genes.use,]

lapply(1:dim(all_dog), function(i){
  make_six_plot(as.character(all_dog[i,1]), as.character(all_dog[i,2]))
})
make_six_plot("CD71", "TFRC")
make_six_plot("CD69", "CD69")
make_six_plot("CD25", "IL2RA")
make_six_plot("CD28", "CD28")
make_six_plot("CD3-1", "CD3E")
make_six_plot("CD279", "PDCD1")

saveRDS(coembed4, file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")

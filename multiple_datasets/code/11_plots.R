library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
source("10a_import_scRNA.R")

# Import main object
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")

# Import matrices
ctrl_scRNA <- import_scRNAseq_cite_stim("ctrl"); dim(ctrl_scRNA)
stim_scRNA <- import_scRNAseq_cite_stim("stim"); dim(stim_scRNA)
ga_control <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_control_ga.rds")
ga_stim <- readRDS("../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_stim_ga.rds")

DimPlot(coembed4, reduction = "umap", group.by = "seurat_clusters", pt.size = .1, split.by = c( 'ca'), label = TRUE, ncol = 2)

DefaultAssay(coembed4) = "ADT"
FeaturePlot(coembed4, reduction = "umap", features = c("CD4-1", "CD19", "CD16", "CD56(NCAM)", "CD8a", "CD14"), pt.size = .3, split.by = c('ca'))

mdf <- data.frame(
  coembed4@reductions$umap@cell.embeddings,
  assay = coembed4@meta.data$orig.ident,
  stim = coembed4@meta.data$stim,
  cluster = coembed4@meta.data$seurat_clusters,
  barcode = substr(rownames(coembed4@meta.data), 1, 16)
)

color_vec <- c("orange2", "orange4", "firebrick", "purple2"); names(color_vec) <-  c("9", "10", "3", "4")
t_cell_vec <- c("#0868AC","#084081", rev(jdb_palette("brewer_marine")[2:7]))
names(t_cell_vec) <- as.character(c(0,1,2,5,6,7,8,11))

set.seed(1)
p1 <- ggplot(shuf(mdf),aes(x=UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.1) + 
  pretty_plot() + 
  facet_wrap( ~ assay +stim , nrow = 1) +
  scale_color_manual(values = c(color_vec, t_cell_vec)) + theme_void() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())
cowplot::ggsave2(p1, file = "../plots/umap_base_viz.png", width = 12, height = 3, dpi = 500)


# import gene <-> protein mapping
gene_mapping <- fread("../data/marker_gene_mapping.tsv")
genes.use <- unique(c(VariableFeatures(ctrl_scRNA), VariableFeatures(stim_scRNA), gene_mapping[[2]]))
length(genes.use)
genes.use[genes.use %ni% rownames(ctrl_scRNA)]
genes.use <- genes.use[genes.use %in% rownames(ctrl_scRNA)]
length(genes.use)

process_pseudo_bulk_changes <- function(clusters, out_name){
  celltype_df <- mdf %>% filter(cluster %in% clusters)
  
  # Quick function for counts per million
  cpm <- function(mat){
    rs <- rowSums(mat)
    return(round(rs/sum(rs)*1000000, 2))
  }
  
  extract_cpm <- function(mat, barcodes){
    (cpm(mat[,substr(colnames(mat), 1, 16) %in% substr(barcodes, 1, 16)]))
  }
  
  # Assemble master table
  atac_control <- extract_cpm(ga_control, celltype_df %>% filter(assay == "ATAC" & stim == "control") %>% pull(barcode))
  atac_stim <- extract_cpm(ga_stim, celltype_df %>% filter(assay == "ATAC" & stim == "stim") %>% pull(barcode))
  
  rna_control <- extract_cpm(ctrl_scRNA@assays$RNA@counts, celltype_df %>% filter(assay == "RNA" & stim == "control") %>% pull(barcode))
  rna_stim <- extract_cpm(stim_scRNA@assays$RNA@counts, celltype_df %>% filter(assay == "RNA" & stim == "stim") %>% pull(barcode))
  
  protein_control <- coembed4@assays$ADT@counts[, rownames(celltype_df)[celltype_df$stim == "control"]] %>% cpm
  protein_stim <- coembed4@assays$ADT@counts[, rownames(celltype_df)[celltype_df$stim == "stim"]] %>% cpm
  
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
  
  cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow = 1), file = paste0("../plots/scatter_DOG_",out_name,".png"), width = 8, height = 5)
  write.table(gene_mapping, file = paste0("../output/DOG_source_",out_name,".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  out_name
} 

process_pseudo_bulk_changes(c(9), "Monocyte")
process_pseudo_bulk_changes(c(4), "Bcell")
process_pseudo_bulk_changes(c(3), "NK")
process_pseudo_bulk_changes(c(0,1,2,5,6,7,8,11), "Tcell")


# Reprocess the t=cell stuff
gene_mapping_t <- fread("../output/DOG_source_Tcell.tsv")
gene_mapping_t_filt <- gene_mapping_t %>% filter(rna_control > 7.5 | rna_stim > 7.5)
gmtf_ss <- gene_mapping_t_filt[,c("log2_atac_change", "log2_rna_change", "log2_protein_change")]; names(gmtf_ss) <- c("ATAC", "RNA", "Protein")
cor(gmtf_ss, use = "pairwise.complete") %>%
  reshape2::melt() -> go_df

p1 <- ggplot(gene_mapping_t_filt,
             aes(x = log2_rna_change, y = log2_protein_change, color = log2_atac_change, label = Marker_name)) +
  geom_point(size = 0.5) + scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot(fontsize = 8) + L_border() +theme(legend.position = "bottom") +
  labs(x = "log2 RNA", y = "log2 protein", color = "log2 atac") +
  scale_x_continuous(limits = c(-5.5, 5.5)) + scale_y_continuous(limits = c(-5.5, 5.5)) +
  theme(legend.position = "none")

cowplot::ggsave2(p1, file = "../plots/dog_scatter_T_for_paper.pdf", 
                 width = 1.45, height = 1.45)

orderf <- c("Protein", "RNA", "ATAC")
go_df$Var1 <- factor(as.character(go_df$Var1), rev(orderf))
go_df$Var2 <- factor(as.character(go_df$Var2), (orderf))

pgrid <- ggplot(go_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + 
  pretty_plot(fontsize = 8) + L_border() +
  theme(legend.position = "none")  +
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0)) +
  labs(x = "", y = "") +
  scale_fill_gradientn(colors = c(jdb_palette("brewer_heat")[c(5:9)]))
cowplot::ggsave2(pgrid, file = "../plots/pgrid_3mode.pdf", 
                 width = 1.85, height = 1.7)

make_six_plot <- function(protein_name, gene_name){
  set.seed(1)
  DefaultAssay(coembed4) = "ADT"
  protein_plot <- FeaturePlot(coembed4, features = protein_name, split.by = "stim", ncol = 1, cells = sample(1:dim(coembed4)[2]),
                              min.cutoff = "q02", max.cutoff = "q98")
  protein_control <- protein_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  protein_stim <- protein_plot[[2]] + theme_void() + theme(legend.position = "none")+ scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  DefaultAssay(coembed4) = "ACTIVITY"
  coembedplot_atac_obj <- coembed4[,coembed4$orig.ident == "ATAC"]
  atac_plot <- FeaturePlot(coembedplot_atac_obj, features = gene_name, split.by = "stim", ncol = 1,cells = sample(1:dim(coembedplot_atac_obj)[2]),
                           min.cutoff = "q02", max.cutoff = "q98")
  atac_control <- atac_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red")) + ggtitle(NULL)
  atac_stim <- atac_plot[[2]] + theme_void() + theme(legend.position = "none")+scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  DefaultAssay(coembed4) = "RNA"
  coembedplot_rna_obj <- coembed4[,coembed4$orig.ident == "RNA"]
  rna_plot <- FeaturePlot(coembedplot_rna_obj, features = gene_name, split.by = "stim", ncol = 1,cells = sample(1:dim(coembedplot_rna_obj)[2]),
                          min.cutoff = "q02", max.cutoff = "q98")
  rna_control <- rna_plot[[1]] + theme_void() + theme(legend.position = "none") + scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  rna_stim <- rna_plot[[2]] + theme_void() + theme(legend.position = "none")+scale_color_gradientn(colors = jdb_palette("ocean_red"))+ ggtitle(NULL)
  
  cowplot::ggsave2(
    cowplot::plot_grid(atac_control, rna_control, protein_control, atac_stim, rna_stim, protein_stim, ncol = 3, scale = 0.9),
    width = 9, height = 6, file = paste0("../plots/dog/", protein_name, "_dog.png"))
  
}

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

all_dog <- gene_mapping[gene_mapping$Gene_symbol %in% genes.use,]
lapply(1:dim(all_dog), function(i){
  make_six_plot(as.character(all_dog[i,1]), as.character(all_dog[i,2]))
})


library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(BuenColors)
library(harmony)

# Pull HQ barcodes from ATAC/RNA embedding
pass_bc_df <- fread("../output/HQ_barcodes_4exps.tsv")
Ac_bcs <- pass_bc_df %>% filter(assay == "ATAC_control") %>% pull(barcode) %>% as.character()
As_bcs <- pass_bc_df %>% filter(assay == "ATAC_stim") %>% pull(barcode) %>% as.character()
Cc_bcs <- pass_bc_df %>% filter(assay == "RNA_control") %>% pull(barcode) %>% as.character()
Cs_bcs <- pass_bc_df %>% filter(assay == "RNA_stim") %>% pull(barcode) %>% as.character()

# Process kite counts
import_kite_counts <- function(library, bcs, tech, bio){
  
  # Import the goodies
  bcs <- substr(bcs, 1, 16)
  mtx <- fread(paste0(library,"featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  
  # Label features and increase stringency to help with embedding
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
asap_control <- import_kite_counts("../data/asapseq/adt/ASAP_ctrl_ADT/",  Ac_bcs, "ASAP", "noStim")
asap_stim <- import_kite_counts("../data/asapseq/adt/ASAP_stim_ADT/",  As_bcs,"ASAP", "Stim")
cite_control <- import_kite_counts("../data/citeseq/adt/ctrl_featurecounts/",  Cc_bcs,"CITE", "noStim")
cite_stim <- import_kite_counts("../data/citeseq/adt/stim_featurecounts/", Cs_bcs,"CITE", "Stim")

pbmc <- merge(merge(merge(asap_control, asap_stim,  add.cell.ids = c("AC", "AS")), cite_control, add.cell.ids= c("", "CC")), cite_stim,add.cell.ids = c("", "CS"))

# Do dimension reduction
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 100, selection.method = "vst")

# Semi-manual curage variable features 
df <- pbmc@assays$ADT@meta.features; df$tag <- rownames(df)
tag_features <- df %>% arrange(desc(vst.variance)) %>% pull(tag) %>% head(150)
pbmc@assays$ADT@var.features <- tag_features

pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunHarmony(pbmc, c("tech", "bio"), assay.use = "ADT")
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = 1:15)

color_vec <- c("orange2", "orange4", "firebrick", "purple2"); names(color_vec) <-  c("9", "10", "3", "4")
t_cell_vec <- c("#0868AC","#084081", rev(jdb_palette("brewer_marine")[2:7]))
names(t_cell_vec) <- as.character(c(0,1,2,5,6,7,8,11))

pass_bc_df$new_assay <- case_when(
  pass_bc_df$assay == "RNA_control" ~ "CITE_noStim",
  pass_bc_df$assay == "RNA_stim" ~ "CITE_Stim",
  pass_bc_df$assay == "ATAC_control" ~ "ASAP_noStim",
  pass_bc_df$assay == "ATAC_stim" ~ "ASAP_Stim"
)
pass_bc_df$barcode_merge <- paste0(pass_bc_df$new_assay, "_",pass_bc_df$barcode)

bc <- paste0(substr(rownames(pbmc@meta.data), nchar(rownames(pbmc@meta.data))-15, nchar(rownames(pbmc@meta.data))), "-1")
seurat_df <- data.frame(pbmc@reductions$umap@cell.embeddings,
                       pbmc@meta.data)
seurat_df$barcode_merge <- paste0(seurat_df$tech, "_", seurat_df$bio, "_", bc)
plot_df <- merge(seurat_df, pass_bc_df, by = "barcode_merge")
plot_df$

set.seed(1)
p1 <- ggplot(shuf(plot_df),aes(x=UMAP_1, y = UMAP_2, color = as.character(cluster))) +
  geom_point(size = 0.1) + 
  pretty_plot() + 
  scale_color_manual(values = c(color_vec, t_cell_vec)) + theme_void() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

p2 <- ggplot(shuf(plot_df),aes(x=UMAP_1, y = UMAP_2, color = bio)) +
  geom_point(size = 0.1) + 
  pretty_plot() + 
  scale_color_manual(values = c("grey", "orange2")) + theme_void() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

p3 <- ggplot(shuf(plot_df),aes(x=UMAP_1, y = UMAP_2, color = tech)) +
  geom_point(size = 0.1) + 
  pretty_plot() + 
  scale_color_manual(values = c("firebrick", "dodgerblue3")) + theme_void() +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text.y = element_blank()) +
  theme(strip.background = element_blank(), strip.text.x = element_blank())

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, nrow = 1, scale = 0.85),
                 file = "../plots/umap_tag_embedding_go.png", width = 12, height = 4, dpi = 500)

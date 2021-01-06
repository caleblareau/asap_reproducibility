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
  colnames(maty) <- paste0(tech, "_", bio, "_", colnames(maty))
  
  # Now create a seurat object
  raw <- CreateSeuratObject(counts = maty,  min.cells = 3, min.features = 10, assay = "ADT"); raw$tech <- tech; raw$bio <- bio
  raw <- NormalizeData(raw, assay = "ADT", normalization.method = "CLR")
  raw <- ScaleData(raw, assay = "ADT")
  
  return(raw)
}

# Import counts 
asap_control <- import_kite_counts("../data/asapseq/dt/ASAP_ctrl_ADT/",  Ac_bcs, "ATAC", "control")
asap_stim <- import_kite_counts("../data/asapseq/adt/ASAP_stim_ADT/",  As_bcs,"ATAC", "stim")
pbmc <- merge(asap_control, asap_stim,  add.cell.ids = c("", ""))

# Do dimension reduction
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, verbose = FALSE)
mat <- pbmc@assays$ADT@scale.data

# A la CiteFuse, train a random forest on the outcome labels to prioritize ADT labels
library(randomForest)
sc_vec <- as.character(pass_bc_df$cluster); names(sc_vec) <- paste0("_", pass_bc_df$assay, "_", substr(pass_bc_df$barcode, 1, 16))
clusters <- sc_vec[colnames(pbmc)]

rf_stim <- randomForest(t(data.matrix(mat)), as.factor(pbmc@meta.data$bio), importance = TRUE)
imp_stim <- importance(rf_stim, type=1, scale = FALSE)

rf_cluster <- randomForest(t(data.matrix(mat)), as.factor(clusters), importance = TRUE)
imp_cluster <- importance(rf_cluster, type=1, scale = FALSE)

importance_df <- data.frame(
  importance_cluster = imp_cluster[,1],
  importance_stim = imp_stim[,1],
  ADT = rownames(mat)
)
importance_df$control <- grepl("sotypeCtrl", importance_df$ADT)

p1 <- ggplot(importance_df %>% arrange(control), aes(x = importance_cluster*100, y = importance_stim*100, color = control, label = ADT)) +
  geom_point() + scale_x_continuous(breaks = c(0,2,4,6)) + scale_y_continuous(breaks = c(0,4,8)) +
  pretty_plot(fontsize = 8) + L_border() + scale_color_manual(values = c("black", "red")) +
  labs(x = "Cluster importance", y = "Stimulation importance") + theme(legend.position = "none")

cowplot::ggsave2(p1, file = "../plots/ADT_importance_ASAP.pdf", width =2, height = 2)

ggplot(importance_df %>% arrange(control), aes(x = importance_cluster*100, y = importance_stim*100, color = control, label = ADT)) +
  geom_text()
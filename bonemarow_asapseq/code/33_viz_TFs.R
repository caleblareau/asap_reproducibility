library(pheatmap)
library(Seurat)
library(matrixStats)
library(BuenColors)

mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Look at the ADT data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
mat <- data.matrix(adtbm@assays$ADT@scale.data)

# Look at the chromVAR data
tfs <- readRDS( "../output/TF_cell_matrix_bagged.rds")
cor_mat <- cor(t(tfs), t(mat), use = "pairwise.complete")

rv  <- rowVars(tfs)
rv_adt  <- rowVars(data.matrix(adt_ss))
names(rv_adt) <- rownames(adt_ss)
threshold <- sort(rv, decreasing = TRUE)[25]
threshold_adt <- sort(rv_adt, decreasing = TRUE)[25]

pp <- pheatmap(cor_mat[(rv >= threshold) & (rowSums(is.na(cor_mat)) == 0),rv_adt>=threshold_adt],
          color = jdb_palette("brewer_spectra"), fontsize = 4)
pdf("../plots/viz_tf_protein.pdf", width = 3.4, height = 3)
pp
dev.off()




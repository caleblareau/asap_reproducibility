library(Seurat)
library(viridis)
library(scales)
library(Matrix)
library(ComplexHeatmap)
library(BuenColors)
library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

# Import ATAC processed data
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Order the cells and the gene / tags
mat <- readRDS("../output/TF_cell_matrix_bagged.rds")
colnames(mat) <- gsub("ASAP_marrow_hg38#", "", colnames(mat))
ordered_cells <- mdf[!is.na(mdf$monocyte_PS),] %>% arrange((monocyte_PS))
tf_mat_ordered <- mat[,ordered_cells$barcode]

# Now group/order/smooth the cell states and make ArchR's nice heatmaps
# see: https://github.com/GreenleafLab/ArchR/blob/5f855d7b7ff3f57eb0d28f312c5ea5d373d27ebd/R/Trajectory.R
n_bins <- 100

groups_tf <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(monocyte_PS >= (idx-1)*multiplier & monocyte_PS < idx*multiplier) %>% pull(barcode)
  rowMeans(tf_mat_ordered[,bcs], na.rm = TRUE)
})

# Filter based on indices
smoothWindow = 11
smooth_groups_tf <- data.matrix((apply((groups_tf), 1, function(x) ArchR:::.centerRollMean((x), k = smoothWindow))))
smooth_groups_minmax_tf <- t(apply(groups_tf, 2, function(x)(x-min(x))/(max(x)-min(x))))


pdf("../plots/monocytic_heatmap_tf.pdf", width = 4, height = 2)
Heatmap(t(smooth_groups_minmax_tf)[c("CEBPB", "GATA1","JUND","SPIB", "KLF1"),],
        col=as.character(jdb_palette("solar_rojos",type="continuous")),
        show_row_names = TRUE, 
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        row_names_gp = gpar(fontsize = 4),
        column_names_gp = gpar(fontsize = 0),
        show_column_names = FALSE)
dev.off()




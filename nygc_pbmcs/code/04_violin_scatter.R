library(data.table)
library(BuenColors)
library(Seurat)
library(dplyr)
source("00_import_kite.R")

sc_LLL <- fread("../data/AvsB/singlecell_LLL.csv.gz") %>% filter(cell_id != "None")
sc_OMNI <- fread("../data/AvsB/singlecell_OMNI.csv.gz")%>% filter(cell_id != "None")

bc_LLL <- sc_LLL$barcode %>% as.character()
bc_OMNI <- sc_OMNI$barcode %>% as.character()

CLR_Seurat <- function(mat){
  pbmc <- CreateSeuratObject(counts = mat, assay = "ADT")
  pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
  pbmc <- ScaleData(pbmc, assay="ADT")
  return(data.matrix(pbmc@assays$ADT@data))
}

# Split into experiment
LLL_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_LLL/featurecounts/", "", FALSE, string_col = "TSA_")
OMNI_kite_A <- import_kite_counts_NYGC("../data/AvsB/TSA_OMNI/featurecounts/", "", FALSE, string_col = "TSA_")
intra_kite_A <- import_kite_counts_NYGC("../data/intra/", "A", TRUE, string_col = "TSA_")

LLL_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_LLL/featurecounts/", "", FALSE, string_col = "TSB_")
OMNI_kite_B <- import_kite_counts_NYGC("../data/AvsB/TSB_OMNI/featurecounts/", "", FALSE, string_col = "TSB_")
intra_kite_B <- import_kite_counts_NYGC("../data/intra/", "B", TRUE, string_col = "TSB_")

LLL_mat <- CLR_Seurat(cbind(LLL_kite_A[bc_LLL,],LLL_kite_B[bc_LLL,]))
OMNI_mat <- CLR_Seurat(cbind(OMNI_kite_A[bc_OMNI,],OMNI_kite_B[bc_OMNI,]))

pal <- colorRampPalette(c("white", "dark green"))
greens_color_ramp <- colorRampPalette(c("white", greens = pal(10)))

plot_marker <- function(mat, marker){
  go_df <- data.frame(
    X = mat[,paste0("TSA_", marker)],
    Y = mat[,paste0("TSB_", marker)]
  )
  plot <- smoothScatter(go_df$X, go_df$Y, nbin = 328, colramp = greens_color_ramp, nrpoints = 0,
                ret.selection = FALSE, xlab="", ylab="", xlim=c(0,4), ylim=c(0,4), labels = FALSE)
  return(plot)
}

markers <- c("CD3", "CD4", "CD8", "CD19", "CD16", "CD56", "CD11c", "CD14", "CD45")
par(mfrow=c(3,3))
par(mar=c(1,1,1,1))
for(i in markers){
  (plot_marker(OMNI_mat, i))
}

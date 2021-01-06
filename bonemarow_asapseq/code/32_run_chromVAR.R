library(Signac)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(chromVARmotifs)
library(dplyr)
library(stringr)
library(SummarizedExperiment)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVARxx)

# Import the HQ cells
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Import peaks x cells and filter for overlap with what's in the gene score / adt pairs
h5_raw <- Read10X_h5("../../../asap_large_data_files/bonemarrow_data/input/asap_marrow_hg38_raw_peak_bc_matrix.h5")
peak_10x_gr <- data.frame(str_split_fixed(rownames(h5_raw), "-|:", 3)) %>% setNames(c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
peak_mat_filt <- h5_raw[,mdf$barcode]

# Set up chromVAR
SE <- SummarizedExperiment(
  rowRanges = peak_10x_gr, 
  assays = list(counts = peak_mat_filt),
  colData = DataFrame(mdf)
)
SE <- addGCBias(SE, BSgenome.Hsapiens.UCSC.hg38)
SE <- filterPeaks(SE, min_fragments = 5)
mm <- matchMotifs(human_pwms_v2, SE, BSgenome.Hsapiens.UCSC.hg38)
dev <- computeDeviations(SE, mm)
bagged <- bagDeviations(dev, cor = 0.7, "human")
dim(dev)
dim(bagged)
tfs <- data.matrix(assays(dev)[["z"]])
saveRDS(tfs, file = "../output/TF_cell_matrix.rds")

tfs2 <- data.matrix(assays(bagged)[["z"]])
rownames(tfs2) <- rowData(bagged)$name
saveRDS(tfs2, file = "../output/TF_cell_matrix_bagged.rds")

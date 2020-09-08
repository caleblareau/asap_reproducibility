library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(Signac)
library(matrixStats)

# Import main objects and annotations
atac_mat <- Read10X_h5("../../pbmc_stimulation_asapseq/data/cellranger/asap_pbmcstim-control_filtered_peak_bc_matrix.h5")
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")
bcdt <- fread('../output/HQ_barcodes_4exps.tsv')
tcellsboo <- c(coembed4@meta.data$seurat_clusters %in% as.character(c(0,1,2,5,6,7,8,11)))
control_cells <- bcdt$barcode[(bcdt$assay %in% c("ATAC_control")) & tcellsboo]
stim_cells <- bcdt$barcode[(bcdt$assay %in% c("ATAC_stim"))& tcellsboo]
stim_cells <- gsub("-1", "-2", stim_cells)

ss_mat_stim <- atac_mat[,colnames(atac_mat) %in% stim_cells]
ss_mat_control <- atac_mat[,colnames(atac_mat) %in% control_cells]
dim(ss_mat_stim)
dim(ss_mat_control)

# Filter for peaks with a CPM > 1
cpm <- function(counts){
  rs <- rowSums(counts)
  cpms <- round(rs/sum(rs) *1000000,1)
  cpms
}
cpms <- cpm(cbind(ss_mat_stim, ss_mat_control))
ss_mat_stim_test <- ss_mat_stim[cpms > 2,]
ss_mat_stim_control <- ss_mat_control[cpms > 2,]

# Function to do differential abundance testing based on permutations
diff_atac_test <- function(mat1, mat2){
  
  mat <- cbind(mat1, mat2)
  dim(mat)
  
  vec <- c(rep(TRUE, dim(mat1)[2]), rep(FALSE, dim(mat2)[2]))
  
  # Permuted values
  permuted <- sapply(1:100, function(i){
    set.seed(i)
    vec_p <- sample(vec)
    p1 <- rowMeans(mat[,vec_p] != 0)
    p2 <- rowMeans(mat[,!vec_p] != 0)
    p1-p2
  })
  
  # Establish baseline characteristics for differential analyses
  diffDF <- data.frame(
    peak_idx = seq(1:dim(mat)[1]),
    p1 = rowMeans(mat[,vec] != 0),
    p2 = rowMeans(mat[,!vec] != 0),
    p = rowMeans(mat != 0),
    mean_permuted = rowMeans(permuted),
    sd_permuted = rowSds(permuted)
  )
  
  # Formalize test statistics
  diffDF %>% 
    mutate(diffP = p1 - p2) %>%
    mutate(Z = ifelse(sd_permuted > 0, (diffP-mean_permuted)/sd_permuted,0)) %>%
    mutate(pvalue = 2*pnorm(-abs(Z)))%>%
    mutate(log10p = -1*log10(pvalue)) %>%
    mutate(FDR = p.adjust(pvalue)) %>% 
    mutate(log10FDR = -1*log10(FDR)) -> diffDF
  
  diffDF
}

diff_atac <- diff_atac_test( ss_mat_stim_control,ss_mat_stim_test)
diff_atac$logFC <- log2((diff_atac$p1 + 0.0001)/(diff_atac$p2+ 0.0001))

saveRDS(diff_atac, file = "../output/DGE_ATAC_Tcell_Activation.rds")

library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(edgeR)

# Import main object
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")
bcdt <- fread('../output/HQ_barcodes_4exps.tsv')
tcellsboo <- c(coembed4@meta.data$seurat_clusters %in% as.character(c(0,1,2,5,6,7,8,11)))
boo_control <- (bcdt$assay %in% c("ATAC_control")) & tcellsboo
boo_stim <- (bcdt$assay %in% c("ATAC_stim"))& tcellsboo
table(boo_stim)
table(boo_control)

# Do protein
control_protein_mat <- t(coembed4@assays$ADT@data)[boo_control,]
stim_protein_mat <- t(coembed4@assays$ADT@data)[boo_stim,]

lapply(1:dim(stim_protein_mat)[2], function(i){
  wt <- wilcox.test(control_protein_mat[,i], stim_protein_mat[,i])
  fc <- log2(mean(stim_protein_mat[,i])/mean(control_protein_mat[,i]))
  data.frame(
    marker = colnames(stim_protein_mat)[i],
    foldchange = fc,
    pvalue = wt$p.value
  )
}) %>% rbindlist() %>% data.frame() -> df
df$padj <- p.adjust(df$pvalue)
table(df$padj < 0.01)
saveRDS(df, file = "../output/DGE_Protein_Tcell_Activation.rds")

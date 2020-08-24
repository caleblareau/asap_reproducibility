library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(edgeR)

# Import main object
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")
bcdt <- fread('../output/HQ_barcodes_4exps.tsv')
tcellsboo <- c(coembed4@meta.data$seurat_clusters %in% as.character(c(0,1,2,5,6,7,8,11)))
boo_control <- (bcdt$assay %in% c("RNA_control")) & tcellsboo
boo_stim <- (bcdt$assay %in% c("RNA_stim"))& tcellsboo
table(boo_stim)
table(boo_control)

# Do RNA first
control_rna_mat <- t(coembed4@assays$RNA@counts)[boo_control,]
stim_rna_mat <- t(coembed4@assays$RNA@counts)[boo_stim,]

run_edgeRQLFdetect_CL <- function(count, condt) {
  
  dge <- DGEList(count, group = condt)
  dge <- calcNormFactors(dge)
  
  # adjust for total expression
  cdr <- scale(colMeans(count > 0))
  cdr2 <- scale(colSums(count))
  
  design <- model.matrix(~ cdr + condt) 
  dge <- estimateDisp(dge, design = design)
  fit <- glmQLFit(dge, design = design)
  qlf <- glmQLFTest(fit)
  tt <- topTags(qlf, n = Inf, sort.by = "none")
  df <- signif(tt@.Data[[1]], 3)
  df$gene <- rownames(df)
  
  # small function to pull CPM
  ps <- function(which_ones){
    rs <- rowSums(count[,which_ones])
    cpms <- round(rs/sum(rs) *1000000,1)[as.character(df$gene)]
    return(cpms)
  } 
  
  # Pull CPM values
  df$control_cpm <- ps(condt == "c")
  df$stim_cpm <- ps(condt == "s")
  
  # Round
  df$logFC <- round(df$logFC,2)
  df$logCPM <- round(df$logCPM,2)
  df
}

condt <- c(
  rep("c", dim(control_rna_mat)[1]),
  rep("s", dim(stim_rna_mat)[1])
)
counts <- t(data.matrix(rbind(control_rna_mat, stim_rna_mat)))
rs <- rowSums(counts)
cpms <- round(rs/sum(rs) *1000000,1)
counts2 <- counts[cpms > 2,]
RNA_DGE <- run_edgeRQLFdetect_CL(counts2, condt)
RNA_DGE_rank <- RNA_DGE %>% arrange(desc(F))
table(RNA_DGE_rank$FDR < 0.01 & (abs(RNA_DGE_rank$logFC) > 0.5))
saveRDS(RNA_DGE_rank, file = "../output/DGE_RNA_Tcell_Activation.rds")

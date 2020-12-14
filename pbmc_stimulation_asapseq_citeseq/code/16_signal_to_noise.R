library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(BuenColors)

# Import main object
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")
bcdt <- fread('../output/HQ_barcodes_4exps.tsv')
tcellsboo <- c(coembed4@meta.data$seurat_clusters %in% as.character(c(0,1,2,5,6,7,8,11)))

# ID cells
boo_control_cite <- (bcdt$assay %in% c("RNA_control")) & tcellsboo
boo_stim_cite <- (bcdt$assay %in% c("RNA_stim"))& tcellsboo
boo_control_asap <- (bcdt$assay %in% c("ATAC_control")) & tcellsboo
boo_stim_asap <- (bcdt$assay %in% c("ATAC_stim"))& tcellsboo

df <- data.frame(
  ag = rownames(coembed4@assays$ADT@counts),
  control = rowSums(coembed4@assays$ADT@counts[,boo_control_asap]),
  stim = rowSums(coembed4@assays$ADT@counts[,boo_stim_asap])
)


ggplot(df, aes(x = control, y = stim, label = ag)) + 
  geom_text() + scale_y_log10() + scale_x_log10()


signal <- c("CD8a", "CD4-1", "CD2", "CD48", "CD45-2", "CD27", "CD47", "CD45-1")
rat <- as.character(df$ag)[grep("sotype", df$ag)]
other <- c("IgD", "CD154", "CD115", "CD40", "IgM", "CD1d", "XCR1", "CD1c")

get_three <- function(boo, assay, condition){
  data.frame(
    signal = colSums(coembed4@assays$ADT@counts[signal,boo]),
    rat = colSums(coembed4@assays$ADT@counts[rat,boo]),
    other = colSums(coembed4@assays$ADT@counts[other,boo]),
    assay, condition
  )
}

combined_df <- rbind(
  get_three(boo_control_cite, "CITE", "Control"),
  get_three(boo_stim_cite, "CITE", "STIM"),
  get_three(boo_control_asap, "ASAP", "Control"),
  get_three(boo_stim_asap, "ASAP", "STIM")
)

ggplot(combined_df, aes(x = condition, y = log10((signal+1)/(other+1)), color = assay)) +
  geom_boxplot() + theme_classic()

ggplot(combined_df, aes(x = condition, y = log10((signal+1)/(rat+1)), color = assay)) +
  geom_boxplot() + theme_classic()

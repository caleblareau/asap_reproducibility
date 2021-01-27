library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(BuenColors)
library(scales)

# Import main object
coembed4 <- readRDS(file = "../../../asap_large_data_files/pbmc_stim_data/output/22July2020_Seurat_Coembed4.rds")
bcdt <- fread('../output/HQ_barcodes_4exps.tsv')

clusters <- unique(sort(bcdt$cluster))

pull_valuesmodality <- function(modality){
  sapply(clusters, function(clustern){
    boo_barcodes <-  (bcdt$assay == modality & bcdt$cluster == clustern)
    rowMeans(coembed4@assays$ADT@data[,boo_barcodes])
  }) %>% reshape2::melt()
}
vals_mat <- data.frame(
  tag = pull_valuesmodality("ATAC_control")[,1],
  cluster = pull_valuesmodality("ATAC_control")[,2],
  ac = pull_valuesmodality("ATAC_control")[,3],
  as = pull_valuesmodality("ATAC_stim")[,3],
  cc = pull_valuesmodality("RNA_control")[,3],
  cs = pull_valuesmodality("RNA_stim")[,3]
)

# Filter for variable proteins
DefaultAssay(coembed4) <- "ADT"
coembed4 <- FindVariableFeatures(coembed4, nfeatures = 50, selection.method = "vst")
tag_df <- coembed4@assays$ADT@meta.features; tag_df$tag <- rownames(tag_df)
tag_features <- tag_df %>% arrange(desc(vst.variance)) %>% pull(tag) %>% head(50)
# Compute correlations
vals_mat[vals_mat[,1] %in% tag_features,c(2:5)] %>%
  cor()

ggplot(vals_mat, aes(x = cs, y = as, label = tag )) +
  geom_text()
max(vals_mat$as)
max(vals_mat$cs)

go_pal = "black"

pdf("../plots/correlation_heat_control.pdf", width = 2.5, height = 3)
smoothScatter(vals_mat$ac,  vals_mat$cc, nbin = 128, colramp = colorRampPalette(c("white", go_pal)),
              nrpoints = 0, ret.selection = FALSE, xlab="", ylab="", xlim = c(0, 3.85),  ylim = c(0, 3.85))
dev.off()
pdf("../plots/correlation_heat_stim.pdf", width = 2.5, height = 3)
smoothScatter(vals_mat$as,  vals_mat$cs, nbin =128, colramp = colorRampPalette(c("white", go_pal)),
              nrpoints = 0, ret.selection = FALSE, xlab="", ylab="", xlim = c(0, 4.7),  ylim = c(0, 4.7))
dev.off()


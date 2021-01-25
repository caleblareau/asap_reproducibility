library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(BuenColors)

pbmc_lll <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")
af <- readRDS("../output/both_LLL_allele_freqs.rds")

consensus_barcodes <- intersect(colnames(af), colnames(pbmc_lll))
af_ss <- af[,consensus_barcodes]
clusters_ordered <- as.character(pbmc_lll@meta.data[consensus_barcodes,"seurat_clusters"])

stimvec <- substr(colnames(af), 18, 18) == 2
stimvec_permuted <- sample(stimvec)
var_meta_df <- data.frame(
  variant = rownames(af),
  mean = rowMeans(af),
  stim_af = (rowMeans(af[,stimvec])),
  control_af = rowMeans(af[,!stimvec])
)

# Do a per-variant anova
var_meta_df$kruskal_pvalue <- sapply(1:dim(var_meta_df)[1], function(i){
  var_df <- data.frame(
    cluster = clusters_ordered,
    variant =af_ss[i,] 
  )
  res.kw <- kruskal.test(variant ~ cluster, data = var_df)
  (res.kw)[[3]]
})
var_meta_df$kruskal_pvalue_adj <- p.adjust(var_meta_df$kruskal_pvalue)
var_meta_df %>% arrange((kruskal_pvalue)) %>% head()

ggplot(var_meta_df, aes(x = -log10(kruskal_pvalue_adj), y = mean*100)) +
  geom_point(size = 0.5) + scale_y_log10(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)*100) +
  labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
  pretty_plot(fontsize = 8) + L_border()

# Visualize specific variants
df <- data.frame(
  pbmc_lll@reductions$wnn.3.umap@cell.embeddings[consensus_barcodes,],
  data.matrix(t(af_ss))
)


pS1 <- ggplot(df %>% arrange((X10761T.C)), aes(x = Uw3_1, y = Uw3_2, color = X10761T.C > 0.1 )) +
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "firebrick")) 

cowplot::ggsave2(pS1, 
                 file = "../plots/mut10761TC_mtDNA_plot.png", height = 1.8, width = 1.8, dpi = 500)

ggplot(df %>% arrange((X10761T.C)), aes(x = Uw3_1, y = Uw3_2, color = X10761T.C > 0.4)) +
  geom_point()
ggplot(df , aes(x = Uw3_1, y = Uw3_2, color = X2299T.C)) +
  geom_point()
ggplot(df , aes(x = Uw3_1, y = Uw3_2, color = X2596G.A >= 0.1)) +
  geom_point()
ggplot(df %>% arrange((X12172A.G)), aes(x = Uw3_1, y = Uw3_2, color = X12172A.G)) +
  geom_point()

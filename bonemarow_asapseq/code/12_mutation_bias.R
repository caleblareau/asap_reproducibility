library(data.table)
library(dplyr)

archr_df <- readRDS(file = "../output/ArchR_main_metadata.rds")
archr_df$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
mSE <- readRDS("../output/mitoMutaitons_cov10_marrow.rds")
archr_cluster_vec <- archr_df$Clusters; names(archr_cluster_vec) <- as.character(archr_df$barcode) 

c_map<- c("Ery", "Ery", rep("Lym", 5), "Prog", # first 8
           rep("Lym", 5), rep("Myeloid", 8))
names(c_map) <- paste0("C", as.character(1:21))
clusters_ordered <- c_map[archr_cluster_vec[as.character(mSE@colData$sample)]]
var_meta_df <- data.frame(rowData(mSE))

# Make easy to interpet supplemental panels
data.frame(
  X13069GA = assays(mSE)[["allele_frequency"]]["13069G>A",],
  X13711GA = assays(mSE)[["allele_frequency"]]["13711G>A",],
  X16260CT = assays(mSE)[["allele_frequency"]]["16260C>T",],
  lineage = clusters_ordered
) %>%
  group_by(lineage) %>% summarize(mX13069GA = mean(X16260CT))

# Do a per-variant anova
var_meta_df$kruskal_pvalue <- sapply(1:dim(var_meta_df)[1], function(i){
  var_df <- data.frame(
    cluster = clusters_ordered,
    variant = assays(mSE)[["allele_frequency"]][i,] 
  )
  res.kw <- kruskal.test(variant ~ cluster, data = var_df)
  (res.kw)[[3]]
})
var_meta_df$kruskal_pvalue_adj <- p.adjust(var_meta_df$kruskal_pvalue)
var_meta_df %>% arrange((kruskal_pvalue)) %>% head()

pMut_mean <- ggplot(var_meta_df, aes(x = -log10(kruskal_pvalue_adj), y = mean*100)) +
  geom_point(size = 0.5) + scale_y_log10(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)*100) +
  labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(pMut_mean, file = "../plots/mutation_bias_mean.pdf", width = 1.8, height= 1.8)


# Pull in UMAP and variants
var_df <- data.frame(barcode = colnames(mSE), clusters_ordered,
                     data.matrix(t(assays(mSE)[["allele_frequency"]])))
plot_df <- left_join(var_df, archr_df, by = "barcode")

pMain1 <- ggplot(shuf(plot_df), aes(x = UMAP1, y = UMAP2, color = X16260C.T > 0.5 )) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "blue")) 

pMain2 <- ggplot(plot_df %>% arrange((X13069G.A)), aes(x = UMAP1, y = UMAP2, color = X13069G.A > 0.05 )) +
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "purple")) 

cowplot::ggsave2(cowplot::plot_grid(pMain1, pMain2, ncol =1), 
                 file = "../plots/two_mito_mutations_MF.png", height = 1.8*2, width = 0.9*2, dpi = 500)

p2 <- ggplot(shuf(plot_df), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = clusters_ordered )) +
  geom_point(size = 0.2)  + pretty_plot() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("purple2", "firebrick", "orange", "forestgreen", "dodgerblue3")) +
  labs(x = "UMAP1", y = "UMAP2")


p5 <- ggplot(plot_df %>% arrange((X13711G.A)), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = X13711G.A > 0.1 )) +
  geom_point(size = 0.2) + pretty_plot()+ theme(legend.position = "bottom") + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "UMAP1", y = "UMAP2")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, p5, ncol = 2), file = "../plots/lineage_bias_boneMarrow_plot.png", width = 8, height = 15)

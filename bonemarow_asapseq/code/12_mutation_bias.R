library(data.table)
library(dplyr)

archr_df <- readRDS(file = "../output/ArchR_main_metadata.rds")
mSE <- readRDS("../output/mitoMutaitons_cov10_marrow.rds")
archr_cluster_vec <- archr_df$Clusters; names(archr_cluster_vec) <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))

c_map<- c("Ery", "Ery", rep("Tcell", 4), 
          rep("Progenitor", 2), rep("Bcell", 4), rep("Myeloid", 8))
names(c_map) <- paste0("C", as.character(1:20))
clusters_ordered <- c_map[archr_cluster_vec[as.character(mSE@colData$sample)]]
var_meta_df <- data.frame(rowData(mSE))

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

p1 <- ggplot(var_meta_df, aes(x = -log10(kruskal_pvalue_adj), y = mean*100)) +
  geom_point() + scale_y_log10(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)*100) +
  labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
  pretty_plot(size = 0.2) + L_border()

# Pull in UMAP
archr_df2 <- data.frame(archr_df, readRDS("../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/archr_proj_analyzed.rds")@embeddings$UMAP$df)
archr_df2$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
var_df <- data.frame(barcode = colnames(mSE), clusters_ordered,
                     t(assays(mSE)[["allele_frequency"]]))
plot_df <- left_join(var_df, archr_df2, by = "barcode")

p2 <- ggplot(shuf(plot_df), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = clusters_ordered )) +
  geom_point(size = 0.2)  + pretty_plot() + theme(legend.position = "bottom") +
  scale_color_manual(values = c("purple2", "firebrick", "orange", "forestgreen", "dodgerblue3")) +
  labs(x = "UMAP1", y = "UMAP2")

p3 <- ggplot(shuf(plot_df), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = X16260C.T > 0.5 )) +
  geom_point(size = 0.2)  + pretty_plot() + theme(legend.position = "bottom") + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "UMAP1", y = "UMAP2")

p4 <- ggplot(plot_df %>% arrange((X13069G.A)), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = X13069G.A > 0.05 )) +
  geom_point(size = 0.2)  + pretty_plot() + theme(legend.position = "bottom") + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "UMAP1", y = "UMAP2")

p5 <- ggplot(plot_df %>% arrange((X13711G.A)), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = X13711G.A > 0.1 )) +
  geom_point(size = 0.2) + pretty_plot()+ theme(legend.position = "bottom") + scale_color_manual(values = c("lightgrey", "firebrick")) +
  labs(x = "UMAP1", y = "UMAP2")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, p3, p4, p5, ncol = 2), file = "../plots/lineage_bias_boneMarrow_plot.png", width = 8, height = 15)

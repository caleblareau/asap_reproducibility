library(data.table)
library(dplyr)

archr_df <- readRDS(file = "../output/ArchR_main_metadata.rds")
mSE <- readRDS("../output/mitoMutaitons_cov10_marrow.rds")
archr_cluster_vec <- archr_df$Clusters; names(archr_cluster_vec) <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
clusters_ordered <- archr_cluster_vec[as.character(mSE@colData$sample)]
var_meta_df <- data.frame(rowData(mSE))

# Do a per-variant anova
var_meta_df$anova_pvalue <- sapply(1:dim(var_meta_df)[1], function(i){
  var_df <- data.frame(
    cluster = clusters_ordered,
    variant = assays(mSE)[["allele_frequency"]][i,]
  )
  res.kw <- kruskal.test(variant ~ cluster, data = var_df)
  (res.kw)[[3]]
})
var_meta_df %>% arrange((anova_pvalue))

ggplot(var_meta_df, aes(x = -log10(anova_pvalue_adj), y = mean*100)) +
  geom_point() + scale_y_log10(breaks = c(0.001, 0.01, 0.1, 0.5)*100) +
  labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
  pretty_plot() + L_border()

# Pull in UMAP
archr_df2 <- data.frame(archr_df, readRDS("../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/archr_proj_analyzed.rds")@embeddings$UMAP$df)
archr_df2$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
var_df <- data.frame(barcode = colnames(mSE), 
                     t(assays(mSE)[["allele_frequency"]]))
plot_df <- left_join(var_df, archr_df2, by = "barcode")

ggplot(plot_df %>% arrange((X204T.C)), aes(x = IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = X204T.C > 0.1 )) +
  geom_point() 

+
  scale_color_gradientn(colors = jdb_palette("brewer_spectra"))

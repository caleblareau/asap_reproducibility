library(Seurat)

mdf <- data.frame(readRDS("../output/ArchR_main_metadata.rds"))
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
adtbm@meta.data$seurat_clusters <- mdf$Clusters
Idents(adtbm) <- mdf$Clusters
FindMarkers(adtbm, ident.1 = "C9") %>% head(10)

mat <- data.matrix(adtbm@assays$ADT@scale.data)

tag_df <- data.frame(
  tag = rownames(mat)
)

tag_df$kruskal_pvalue <- sapply(1:dim(mat)[1], function(i){
  var_df <- data.frame(
    cluster = mdf$Clusters,
    tag_value = mat[i,]
  )
  res.kw <- kruskal.test(tag_value ~ cluster, data = var_df)
  res.kw[[1]]
})
tag_df %>% arrange(desc(kruskal_pvalue))

dfd <- data.frame(
  cluster = mdf$Clusters,
  t(mat)
) %>% group_by(cluster) %>% summarize_all(mean)


dfd[,c("cluster", "CD4.1","CD34","CD8",  "CD71", "CLEC12A", "CD31", "CD56.NCAM.","CD3.1", "CD127")]

umap_df <- mdf[,c("UMAP1", "UMAP2")]
plot_df <- data.frame(umap_df, t(mat))

plot_factor <- function(v){
  plot_df_one <- data.frame(
    umap_df,
    value = mat[v,]
  )
  p1 <- ggplot(plot_df_one %>% arrange((value)), aes(x = UMAP1, y = UMAP2, color = value)) +
    geom_point(size = 0.2) + ggtitle(v)+ pretty_plot() + L_border() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("samba_color"), limits = c(0,3), oob = scales::squish)
  v2 <- gsub("/", "_", v)
  cowplot::ggsave2(p1, file = paste0("../plots/full_umap/", v2, ".pdf"), width = 4, height = 4.2)
}
lapply(rownames(mat), plot_factor)

library(Seurat)

mdf <- data.frame(readRDS("../output/ArchR_main_metadata.rds"))
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
adtbm@meta.data$seurat_clusters <- mdf$Clusters
Idents(adtbm) <- mdf$Clusters

mat <- data.matrix(adtbm@assays$ADT@scale.data)

umap_df <- mdf[,c("UMAP1", "UMAP2")]
umap_mat <- data.matrix(umap_df)
rownames(umap_mat) <- colnames(adtbm)
adtbm[["umaparchr"]] <- CreateDimReducObject(embeddings = umap_mat, key = "AU_", assay = DefaultAssay(adtbm))

plot_df <- data.frame(umap_df, t(mat))

v <- "CD4-1"
plot_df_one <- data.frame(
  umap_df,
  value = mat[v,]
)  %>% arrange(value) %>% ggplot(aes(x = UMAP1, y = UMAP2, color = value)) +
  geom_point(size = 0.2) + theme_void() + 
  scale_color_gradientn(colors = jdb_palette("samba_color"), limits = c(0,3), oob = scales::squish) -> p_p_for_legend
#cowplot::ggsave2(g_legend(p_p_for_legend), file = "../plots/samba_color_legend.pdf", width = 1, height = 2)

# Function to make little squares for big panel
plot_factor_shotgun <- function(v){
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
#lapply(rownames(mat), plot_factor_shotgun)

# Make selected plot for main figure
plot_factor_qs_m <- function(v){
  plot_df_one <- data.frame(
    umap_df,
    value = mat[v,]
  )
  qs <- quantile(plot_df_one$value, c(0.01, 0.99))
  q1 <- qs[1]
  q99 <- qs[2]
  
  p1 <- ggplot(plot_df_one %>% arrange((value)), aes(x = UMAP1, y = UMAP2, color = value)) +
    geom_point(size = 0.2) + theme_void() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("solar_extra"), limits = c(q1,q99), oob = scales::squish)
  return(p1)
}

plot_factor_qs <- function(v){
   p1 <- FeaturePlot(adtbm, features = c(v), pt.size = 0.2, reduction = "umaparchr", min.cutoff = "q2", max.cutoff = "q98")+ 
     theme_void() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("solar_extra"))
  return(p1)
}

cowplot::ggsave2(
  cowplot::plot_grid(
    plot_factor_qs_m("CLEC12A"),
    plot_factor_qs_m("CD31"),
    plot_factor_qs_m("CD2"),
    plot_factor_qs_m("CD9"),
    plot_factor_qs_m("CD38-1"),
    plot_factor_qs_m("CD56(NCAM)"), ncol = 3, scale = 0.8
  ), file = "../plots/adt_small_panels_supplement_q1_q99.png", width = 2.25*4, height = 1.5*4, dpi = 300)

c("CLEC12A", "CD31", "CD2", "CD9", "CD38-1", "CD56(NCAM)") %in% rownames(mat)

cowplot::ggsave2(
  cowplot::plot_grid(
    plot_factor_qs_m("CD34"),
    plot_factor_qs_m("CD71"),
    plot_factor_qs_m("CD19"),
    plot_factor_qs_m("CD14"),
    plot_factor_qs_m("CD3-2"),
    plot_factor_qs_m("CD4-1"), ncol = 3, scale = 0.8
  ), file = "../plots/adt_small_panels_sa_q1_q99.png", width = 2.25*4, height = 1.5*4, dpi = 300)


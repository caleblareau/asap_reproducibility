library(data.table)
library(BuenColors)
library(Seurat)

mdf <- data.frame(readRDS("../output/ArchR_main_metadata.rds"))
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
adtbm@meta.data$seurat_clusters <- mdf$Clusters
Idents(adtbm) <- mdf$Clusters

# Handle totalVI
dat <- fread("../output/13june2020_citeseq_totalVI_output.csv.gz")
dat$batch <- substr(dat$barcode, 18,18)
dat$total_tags <- rowSums(data.matrix(dat[,seq(28, 472, 2)-1]))

mdf <- merge(dat[,c(1,2,3,4,25,26)], data.frame(barcode = colnames(adtbm), t(adtbm@assays$ADT@scale.data)), by = "barcode")

ggplot(mdf, aes(x = UMAP1, y = UMAP2, color = as.factor(leiden_totalVI))) + 
  geom_point() + pretty_plot() + theme(legend.position = "none") +
  scale_color_manual(values = jdb_palette("corona"))


plot_factor <- function(v){
  plot_df_one <- data.frame(
    mdf[,c("UMAP1", "UMAP2")],
    value = mdf[[v]]
  )
  p1 <- ggplot(plot_df_one %>% arrange((value)), aes(x = UMAP1, y = UMAP2, color = value)) +
    geom_point(size = 0.2) + ggtitle(v)+ pretty_plot() + L_border() + theme(legend.position = "none") +
    scale_color_gradientn(colors = jdb_palette("samba_color"), limits = c(0,3), oob = scales::squish)
  v2 <- gsub("/", "_", v)
  cowplot::ggsave2(p1, file = paste0("../plots/total_vi_umap/", v2, ".pdf"), width = 4, height = 4.2)
}
lapply(colnames(mdf)[5:246], plot_factor)



library(BuenColors)
library(dplyr)
library(ggrastr)

# Set up color annotation / reannotion convention
cluster_reanno <- c("E1", "E2", "T1", "T2", "T3", "T4", "T5",
                    "Proj", "pDC", "B1", "B2", "B3", "B4", "Baso",
                    "M1", "M2", "M3", "M4", "M5", "M6", "M7")
names(cluster_reanno) <- paste0("C", as.character(1:21))

color_vec <- c("#FB6A4A", "#A50F15",  jdb_palette("brewer_marine")[c(3,5,6,7,8)],
               "green4", "#094078", "#581F99", "#6D49AB", "#8880C4", "#9DA8D7", "#CD8500",
               jdb_palette("brewer_fire")[c(4,2,5,6,7,8,9)])
names(color_vec) <- unname(cluster_reanno)

mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))
mdf$nc <- cluster_reanno[as.character(mdf$Clusters)]

set.seed(4)
pEmbed <- ggplot(shuf(mdf), aes(x = UMAP1, y = UMAP2, color = nc)) +
  geom_point_rast(raster.dpi = 1000, size = 0.1) + theme_void(base_size = 6) + labs(color = "Cluster") +
  scale_color_manual(guide = guide_legend(ncol = 3), values = color_vec)
cowplot::ggsave2(pEmbed, file = "../plots/panelB_embedding_basic.pdf", width = 3.7, height = 2)

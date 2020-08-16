library(data.table)
library(BuenColors)
library(dplyr)

dat <- fread("../output/totalVI_output_bonemarrow_asap.csv.gz")

# Import ATAC processed data
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

two_clusters <- merge(mdf[,c("barcode", "nc", "UMAP1", "UMAP2")], dat[,c("barcode","leiden_totalVI", "UMAP1", "UMAP2")], by = "barcode")

set.seed(1)
p_VIumap_VIclusters <- ggplot(shuf(two_clusters), aes(x = UMAP1.y, y = UMAP2.y, color = as.factor(leiden_totalVI))) +
  geom_point(size = 0.1) + theme_void() + labs(color = "Cluster") +
  scale_color_manual(values = jdb_palette("corona")) + theme(legend.position = "none")
p_VIumap_ArchRclusters <- ggplot(shuf(two_clusters), aes(x = UMAP1.y, y = UMAP2.y, color = nc)) +
  geom_point(size = 0.1) + theme_void() + labs(color = "Cluster") +
  scale_color_manual(values = color_vec) + theme(legend.position = "none")
p_ArchrRumap_VIclusters <- ggplot(shuf(two_clusters), aes(x = UMAP1.x, y = UMAP2.x, color = as.factor(leiden_totalVI))) +
  geom_point(size = 0.1) + theme_void() + labs(color = "Cluster") +
  scale_color_manual(values = jdb_palette("corona")) + theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(p_VIumap_VIclusters, p_VIumap_ArchRclusters,p_ArchrRumap_VIclusters, nrow = 1, scale = 0.9),
                 file = "../plots/compare_archr_totalVI_embeddings.pdf", width = 5.4, height = 1.8, dpi = 500)



#--- look at confusion matrices

confusion_df2 <- two_clusters %>% 
  group_by(leiden_totalVI, nc) %>% summarize(total = n())
confusion_df2 <- confusion_df2 %>% group_by(nc) %>% mutate(fill_val = total / sum(total))
confusion_df2$nc <- factor(as.character(confusion_df2$nc),
                           levels = sort(unname(cluster_reanno), decreasing = TRUE))

# Order and plot
confusion_df3 <- confusion_df2 %>% arrange(desc(nc)) %>% arrange(desc(fill_val > 0.1)) %>% arrange(desc(fill_val > 0.3))
confusion_df3$leiden <- factor(as.character(confusion_df3$leiden_totalVI), 
                               levels = unique(as.character(confusion_df3$leiden_totalVI)))
ggplot(confusion_df3, aes(x = leiden, y = nc, label = total, fill = fill_val*100)) +
  geom_tile(color = NA) +
  pretty_plot(fontsize = 8) + scale_fill_gradientn(colors = jdb_palette("brewer_heat")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  scale_color_manual(guide = "none", values = c("black", "white")) + 
  labs(x = "totalVI clusters", y = "ArchR clusters", fill = "% of ArchR cluster")

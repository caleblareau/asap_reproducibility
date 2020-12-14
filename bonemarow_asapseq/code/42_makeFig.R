library(BuenColors)

mdf <- readRDS("../output/Seurat_Proteinprojections_allMeta.rds")
dd <- mdf %>% group_by(predictions_l2, archr_cluster) %>%
  summarize(count = n()) %>% ungroup() %>% group_by(predictions_l2) %>% mutate(prop = count/sum(count)*100)

p1 <- ggplot(dd, aes(x = archr_cluster, y = predictions_l2, fill =  prop)) +
  geom_tile() + pretty_plot() + L_border() +
  scale_fill_gradientn(colours = jdb_palette("brewer_heat")) +
  theme(legend.position = "bottom") +
  scale_y_discrete(limits = rev(unique(dd$predictions_l2)))+ 
  labs(x = "ArchR Cluster", y = "Seurat Projection", fill= "% of Seurat projection")

d2 <- mdf %>% group_by(predictions_l2, archr_cluster) %>%
  summarize(count = n()) %>% ungroup() %>% group_by(archr_cluster) %>% mutate(prop = count/sum(count)*100)

p2 <- ggplot(d2, aes(x = archr_cluster, y = predictions_l2, fill =  prop)) +
  geom_tile() + pretty_plot() + L_border() +
  scale_fill_gradientn(colours = jdb_palette("brewer_heat")) +
  theme(legend.position = "bottom") +
  scale_y_discrete(limits = rev(unique(d2$predictions_l2)))+ 
  labs(x = "ArchR Cluster", y = "Seurat Projection", fill= "% of ArchR cluster")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow =1), file = "../plots/heatmaps_proportions.png",
                 width = 12, height = 5)
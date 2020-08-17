library(ArchR)
doubScores <- addDoubletScores(
  input = readRDS("../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/archr_proj_analyzed.rds")
)
df <- data.frame(DS = doubScores@cellColData$DoubletScore, 
                 DE = doubScores@cellColData$DoubletEnrichment, 
                 readRDS("../output/ArchR_main_metadata.rds"))

df$DE2 <- ifelse(df$DE > 3, 3, df$DE)
p1 <- ggplot(shuf(df) , aes(x = UMAP1, y = UMAP2, color = DE2)) +
  geom_point(size = 0.3) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, 
                 file = "../plots/umap_doublet.png", width = 1.8*4, height = 1.8*4, dpi = 300)

p2 <- qplot(df$DE, bins = 50) +
  labs(x = "Doublet Enrichment", y = "Count") +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p2, file = "../plots/doublet_histogram.pdf", width = 1.8, height = 1.8)

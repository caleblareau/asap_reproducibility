library(data.table)
library(BuenColors)
library(dplyr)

dat <- fread("../output/totalVI_output_bonemarrow_asap.csv.gz")

# Import ATAC processed data
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

two_clusters <- merge(mdf[,c("barcode", "Clusters")], dat[,c("barcode","leiden_totalVI")])

confusion_df2 <- two_clusters %>% 
  group_by(leiden_totalVI, Clusters) %>% summarize(total = n())
confusion_df2 <- confusion_df2 %>% group_by(Clusters) %>% mutate(fill_val = total / sum(total))

ggplot(confusion_df2, aes(x = as.factor(leiden_totalVI), y = Clusters, label = total, fill = fill_val*100, color = total > 1000)) +
  geom_tile(color = "black") +
  pretty_plot(fontsize = 8) + scale_fill_gradientn(colors = jdb_palette("brewer_heat")) +
  theme(legend.position = "bottom") + 
  scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) +
  scale_color_manual(guide = "none", values = c("black", "white")) + 
  labs(x = "totalVI clusters", y = "ArchR clusters", fill = "% of ArchR cluster")

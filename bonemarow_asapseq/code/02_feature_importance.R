library(Seurat)
library(viridis)
library(scales)

# Import ATAC processed data
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Look at the ADT data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
mat <- data.matrix(adtbm@assays$ADT@scale.data)

# A la CiteFuse, train a random forest on the outcome labels to prioritize ADT labels
library(randomForest)

rf_bm <- randomForest(t(data.matrix(mat)), as.factor(mdf$Clusters ), importance = TRUE)
imp_bm <- importance(rf_bm, type=1, scale = FALSE)

importance_df <- data.frame(
  importance = imp_bm[,1],
  ADT = rownames(mat)
) %>% arrange(desc(importance)) %>% mutate(rank = 1:n())
importance_df$control <- grepl("sotypeCtrl", importance_df$ADT)

p_importance_adt <- ggplot(importance_df %>% arrange(control), aes(x = rank, y = importance*100, color = control)) +
  geom_point(size = 0.5) +  scale_y_continuous(breaks = c(0,2.5, 5)) + 
  pretty_plot(fontsize = 8) + L_border() + scale_color_manual(values = c("black", "red")) +
  labs(x = "Ranked ADTs", y = "ADT importance") + theme(legend.position = "none")

cowplot::ggsave2(p_importance_adt, file = "../plots/ADT_importance_bonemarrow.pdf", width = 1.8, height = 1.8)


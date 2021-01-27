library(SummarizedExperiment)
library(Matrix)

# Function that quickly computes the allele frequency matrix from a summarized experiment mgatk object
computeAFMutMatrix <- function(SE){
  cov <- assays(SE)[["coverage"]]+ 0.001
  ref_allele <- as.character(rowRanges(SE)$refAllele)
  
  getMutMatrix <- function(letter){
    mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
    rownames(mat) <- paste0(as.character(1:dim(mat)[1]), toupper(ref_allele), ">", letter)
    return(mat[toupper(ref_allele) != letter,])
  }
  
  rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
  
}

dir <- "../../../asap_large_data_files/broad_experiment_pbmcs/mgatk_input/"
af_A <- t(computeAFMutMatrix(readRDS(paste0(dir, "ASAP_A_v12_hg38-mtMask_mgatk.rds")))[c("5979G>A", "12763G>A", "4799C>T"),])
af_B <- t(computeAFMutMatrix(readRDS(paste0(dir, "ASAP_B_v12_hg38-mtMask_mgatk.rds")))[c("5979G>A", "12763G>A", "4799C>T"),])
af_C <- t(computeAFMutMatrix(readRDS(paste0(dir, "ASAP_C_v12_hg38-mtMask_mgatk.rds")))[c("5979G>A", "12763G>A", "4799C>T"),])

mito_df <- data.frame(
  Barcode = c(rownames(af_A), gsub("-1", "-2", rownames(af_B)), gsub("-1", "-3", rownames(af_C))),
  data.matrix(rbind(af_A, af_B, af_C)),
  experiment = c(rep("A", dim(af_A)[1]), rep("B", dim(af_B)[1]), rep("C", dim(af_C)[1]))
)

cr_df <- data.frame(
  fread("../../../asap_large_data_files/broad_experiment_pbmcs/cellranger/analysis/tsne/2_components/projection.csv"),
  fread("../../../asap_large_data_files/broad_experiment_pbmcs/cellranger/analysis/clustering/graphclust/clusters.csv")
)

mdf <- merge(cr_df, mito_df, by = "Barcode")
saveRDS(mdf, file = "../output/BroadExperiment_cellranger_ouput_withMito.rds")

library(cowplot)
p1m1 <- ggplot(mdf %>% filter(experiment == "A") %>% arrange((X5979G.A)), aes(x = TSNE.1, y = TSNE.2, color = X5979G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p1m2 <- ggplot(mdf %>% filter(experiment == "B") %>% arrange((X5979G.A)), aes(x = TSNE.1, y = TSNE.2, color = X5979G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p1m3 <- ggplot(mdf %>% filter(experiment == "C") %>% arrange((X5979G.A)), aes(x = TSNE.1, y = TSNE.2, color = X5979G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p2m1 <- ggplot(mdf %>% filter(experiment == "A") %>% arrange((X4799C.T)), aes(x = TSNE.1, y = TSNE.2, color = X4799C.T)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p2m2 <- ggplot(mdf %>% filter(experiment == "B") %>% arrange((X4799C.T)), aes(x = TSNE.1, y = TSNE.2, color = X4799C.T)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p2m3 <- ggplot(mdf %>% filter(experiment == "C") %>% arrange((X4799C.T)), aes(x = TSNE.1, y = TSNE.2, color = X4799C.T)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p3m1 <- ggplot(mdf %>% filter(experiment == "A") %>% arrange((X12763G.A)), aes(x = TSNE.1, y = TSNE.2, color = X12763G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p3m2 <- ggplot(mdf %>% filter(experiment == "B") %>% arrange((X12763G.A)), aes(x = TSNE.1, y = TSNE.2, color = X12763G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))

p3m3 <- ggplot(mdf %>% filter(experiment == "C") %>% arrange((X12763G.A)), aes(x = TSNE.1, y = TSNE.2, color = X12763G.A)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_gradientn(colors = c("lightgrey", "firebrick"))


ggsave(cowplot::plot_grid(p1m1, p2m1, p1m2, p2m2,  p1m3, p2m3, ncol = 2), 
       file = "plots/p3_m3.png", width = 3.66, height = 5, dpi = 500)

mdf$Cluster <- as.character(mdf$Cluster)
pc1 <- ggplot(mdf %>% filter(experiment == "A"), aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_manual(values = jdb_palette("corona"))
pc2 <- ggplot(mdf %>% filter(experiment == "B"), aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_manual(values = jdb_palette("corona"))
pc3 <- ggplot(mdf %>% filter(experiment == "C"), aes(x = TSNE.1, y = TSNE.2, color = Cluster)) +
  geom_point(size = 0.05) + theme_nothing()  + theme(legend.position = "none") +
  scale_color_manual(values = jdb_palette("corona"))
ggsave(cowplot::plot_grid(pc1, pc2, pc3, ncol = 1), 
       file = "plots/cell_clusters_scATAC.png", width = 1.666, height = 5, dpi = 500)


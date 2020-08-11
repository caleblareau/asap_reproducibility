library(ArchR)
library(dplyr)

# Establish ArchR project
proj <- ArchRProject(
  ArrowFiles = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38.arrow", 
  outputDirectory = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38",
  copyArrows = FALSE
)

# Filter cells
all_cells <- proj$cellNames
valid_cells <- paste0("ASAP_marrow_hg38#", fread("../data/barcodes/step3_ADThq.tsv", header = FALSE)[[1]])

subset_cells <- (all_cells %in% valid_cells) 
keep_cells <- all_cells[subset_cells]
proj <- subsetCells(ArchRProj = proj, cellNames = keep_cells)

# Do ArchR things
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addImputeWeights(proj)
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", minDist = 0.4, force = TRUE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI", resolution = 1, force = TRUE)

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

# Pull out UMAP for manual analysis
umap_df <- data.frame(proj@embeddings$UMAP@listData$df); colnames(umap_df) <- c("UMAP1", "UMAP2")
meta_df <- data.frame(
  proj@cellColData,
  umap_df
)
ggplot(meta_df %>% dplyr::filter(Clusters == "C8" & UMAP1 < 0), aes(x = UMAP1, y = UMAP2)) + geom_point() +
  geom_vline(xintercept = -4, color = "firebrick") +
  geom_vline(xintercept = -4.7, color = "dodgerblue3") 
  

proj@cellColData$Clusters_for_ery_PS <- case_when(
  meta_df$Clusters == "C1" ~ "C1",
  meta_df$Clusters == "C2" ~ "C2",
  meta_df$Clusters == "C8" & meta_df$UMAP1 > -4.7 ~ "C0e",
  TRUE ~ "other"
)

proj@cellColData$Clusters_for_mono_PS <- case_when(
  meta_df$Clusters == "C20" ~ "C20",
  meta_df$Clusters == "C21" ~ "C21",
  meta_df$Clusters == "C8" & meta_df$UMAP1 < -4 ~ "C0m",
  TRUE ~ "other"
)

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters_for_mono_PS", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters_for_ery_PS", embedding = "UMAP")

# Add a pseudotime
proj <- addTrajectory(
  ArchRProj = proj, 
  name = "monocyte_PS", 
  groupBy = "Clusters_for_mono_PS",
  trajectory = c( "C0m", "C20", "C21"), 
  embedding = "UMAP", 
  force = TRUE
)

proj <- addTrajectory(
  ArchRProj = proj, 
  name = "erythroid_PS", 
  groupBy = "Clusters_for_ery_PS",
  trajectory = c( "C0e", "C1", "C2"), 
  embedding = "UMAP", 
  force = TRUE
)
library(BuenColors)
p_mps <- plotTrajectory(proj, trajectory = "monocyte_PS", colorBy = "cellColData", name = "monocyte_PS")
cowplot::ggsave2(p_mps[[1]] + theme_void() + theme(legend.position = "none") + ggtitle(""), file = "../plots/pseudotime_monocyte.pdf",
                 width = 6, height = 6)
p_eps <- plotTrajectory(proj, trajectory = "erythroid_PS", colorBy = "cellColData", name = "erythroid_PS")
cowplot::ggsave2(p_eps[[1]] + theme_void() + theme(legend.position = "none") + ggtitle(""), file = "../plots/pseudotime_erythroid.pdf",
                 width = 6, height = 6)

meta_df <- data.frame(
  proj@cellColData,
  umap_df
)
saveRDS(meta_df, file = "../output/ArchR_main_metadata.rds")
saveRDS(proj, file = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/archr_proj_analyzed.rds")

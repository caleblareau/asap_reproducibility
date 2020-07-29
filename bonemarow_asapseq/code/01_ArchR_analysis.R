library(ArchR)

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

# Visualize marker genes
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1",  #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A", "CD4" #TCells
)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

#Rearrange for grid plotting
p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# Add a pseudotime
proj <- addTrajectory(
  ArchRProj = proj, 
  name = "myeloidPS", 
  groupBy = "Clusters",
  trajectory = c( "C8", "C20", "C21"), 
  embedding = "UMAP", 
  force = TRUE
)
p <- plotTrajectory(proj, trajectory = "myeloidPS", colorBy = "cellColData", name = "myeloidPS")
p[[1]]

umap_df <- data.frame(proj@embeddings$UMAP@listData$df); colnames(umap_df) <- c("UMAP1", "UMAP2")
meta_df <- data.frame(
  proj@cellColData,
  umap_df
)
saveRDS(meta_df, file = "../output/ArchR_main_metadata.rds")
saveRDS(proj, file = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/archr_proj_analyzed.rds")

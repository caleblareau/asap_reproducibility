library(ArchR)

path <- "../../../asap_large_data_files/pbmc_stim_data/input/"
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = paste0(path, c("stim_fragments.tsv.gz","control_fragments.tsv.gz")),
  sampleNames = c("Stim", "Control"),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Establish ArchR project
path_arrow <- "../../../asap_large_data_files/pbmc_stim_data/output/archr_pbmc_stim/"
proj <- ArchRProject(
  ArrowFiles = paste0(path_arrow, c("Control.arrow", "Stim.arrow")), 
  outputDirectory = paste0(path_arrow, "ASAP_PBMCstim_hg38"),
  copyArrows = FALSE
)

proj <- addDoubletScores(input = proj,
                         k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                         knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                         LSIMethod = 1
)
proj <- filterDoublets(ArchRProj = proj, cutScore = 30)

# Subset to cells with valid ADT data
all_cells <- proj$cellNames
valid_adt_cells <- c(
  paste0("Control#", colnames(readRDS("../output/adt_mat/ASAP_ctrl.rds")), "-1"),
  paste0("Stim#", colnames(readRDS("../output/adt_mat/ASAP_stim.rds")), "-1")
)

# Filter cells in ArchR project
subset_cells <- (all_cells %in% valid_adt_cells) & (proj@cellColData$DoubletScore < 30)
table(subset_cells)
keep_cells <- all_cells[subset_cells]
proj <- subsetCells(ArchRProj = proj, cellNames = keep_cells)

# Compute cell-cell similarity using MAGIC
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", force = TRUE)
proj <- addImputeWeights(proj)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony", name = "UMAP", force = TRUE,
                nNeighbors = 20, minDist = 0.2)

proj <- addClusters(input = proj, reducedDims = "Harmony", name = "Clusters",
                    force = TRUE, resolution = 0.4)

plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAP")

# Visualize marker genes
markerGenes  <- c(
  "MS4A1", "CD19", "PTPRC", #B-Cell etc
  "CD14", "MPO", "FCGR3A",  #Monocytes
  "CD3D", "CD8A", "CD4", #TCells
  "IFNG", "CD69", "NCAM1"
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
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
saveRDS(proj, file = paste0(path_arrow, "ArchR_PBMCs_stim.rds"))



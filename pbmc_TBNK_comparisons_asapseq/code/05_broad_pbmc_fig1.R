#### Figure 1 #####
# Code from Eleni
### making the ArchR object for the broad experiment

input.file.list = "../../../asap_large_data_files/broad_experiment_pbmcs/cellranger_individual/ASAP_B_v12_hg38_fragments.tsv.gz"
sample.name.ls = "broad"

ArrowFiles <- createArrowFiles(
  inputFiles = input.file.list,
  sampleNames = sample.name.ls,
  filterTSS = 4, 
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = ".",
  copyArrows = TRUE
)

### importing scaled ADT values from Broad's seurat object

broad.so = readRDS("../output/Broad_ASAPB_processed_ADT_object.rds")

ADT=broad.so@assays$ADT@scale.data

## fixing the _2 in barcode names
tmp = strsplit(colnames(ADT), "-")
tmp2 = sapply(tmp, function(dat) substr(dat[1],1,16))
tmp3 <- paste(tmp2, "1", sep = "-")
colnames(ADT) = tmp3

colnames(ADT) <- paste0("broad#", colnames(ADT))
cells.to.keep = intersect(proj$cellNames, colnames(ADT))
subset <- subsetCells(ArchRProj = proj, cellNames = cells.to.keep)

adt.matched = ADT[, match(subset$cellNames, colnames(ADT))]

subset <- addIterativeLSI(ArchRProj = subset, useMatrix = "TileMatrix", name = "IterativeLSI")
subset <- addImputeWeights(subset)
subset <- addUMAP(ArchRProj = subset, reducedDims = "IterativeLSI", minDist = 0.4, force = TRUE)
subset <- addClusters(input = subset, reducedDims = "IterativeLSI")

subset@cellColData$CD45 = adt.matched[9,]
subset@cellColData$CD8 = adt.matched[8,]
subset@cellColData$CD14 = adt.matched[7,]
subset@cellColData$CD56 = adt.matched[6,]
subset@cellColData$CD11c = adt.matched[5,]
subset@cellColData$CD4 = adt.matched[4,]
subset@cellColData$CD16 = adt.matched[3,]
subset@cellColData$CD3 = adt.matched[2,]
subset@cellColData$CD19 = adt.matched[1,]

### importing heteroplasmy for 3 variants

mito = readRDS("../output/BroadExperiment_cellranger_ouput_withMito.rds")
b = mito[which(mito$experiment=="B"),]
tmp = strsplit(b$Barcode, "-")
tmp2 = sapply(tmp, function(dat) substr(dat[1],1,16))
tmp3 <- paste(tmp2, "1", sep = "-")
b$Barcode.1 = tmp3
mito = t(mito)
mito = b[,c(4,6,7,8)]
rownames(mito) = mito$Barcode.1
mito$Barcode.1 = NULL
temp = t(mito)
colnames(temp) <- paste0("broad#", colnames(temp))
mito.matched = temp[, match(subset$cellNames, colnames(temp))]


subset@cellColData$X5979G.A = mito.matched[1,]
subset@cellColData$X12763G.A = mito.matched[2,]
subset@cellColData$X4799C.T = mito.matched[3,]

x59bc = subset$cellNames[which(subset@cellColData$X5979G.A>0)]
x47bc = subset$cellNames[which(subset@cellColData$X4799C.T>0)]

# panel e

p1 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD4", embedding = "UMAP", quantCut = c(0.1,0.9))
p2 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD8", embedding = "UMAP", quantCut = c(0.1,0.9))
p5 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD14", embedding = "UMAP", quantCut = c(0.1,0.9))
p6 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD16", embedding = "UMAP", quantCut = c(0.1,0.9))
plot_grid(p1, p2, p6, p5, nrow = 2)

# panel f

plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "X5979G.A", embedding = "UMAP",
              highlightCells = x59bc, plotAs = "X5979G.A",
              pal=ArchRPalettes$greyMagma, rastr = FALSE, size = 0.01)
plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "X4799C.T", embedding = "UMAP",
              highlightCells = x47bc, plotAs = "X4799C.T",
              pal=ArchRPalettes$greyMagma, rastr = FALSE, size = 0.01)


#### Suppl Fig 1 ####

# panel h

plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

# panel i

p1 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD4", embedding = "UMAP", quantCut = c(0.1,0.9))
p2 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD8", embedding = "UMAP", quantCut = c(0.1,0.9))
p3 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD19", embedding = "UMAP", quantCut = c(0.1,0.9))
p4 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD11c", embedding = "UMAP", quantCut = c(0.1,0.9))
p5 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD14", embedding = "UMAP", quantCut = c(0.1,0.9))
p6 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD16", embedding = "UMAP", quantCut = c(0.1,0.9))
p7 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD45", embedding = "UMAP", quantCut = c(0.1,0.9))
p8 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD3", embedding = "UMAP", quantCut = c(0.1,0.9))
p9 = plotEmbedding(ArchRProj = subset, colorBy = "cellColData", name = "CD56", embedding = "UMAP", quantCut = c(0.1,0.9))
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)

# panel j

markerGenes  <- c(
  "CD4", 
  "CD8A", 
  "CD19",
  "ITGAX",
  "CD14",
  "FCGR3A",
  "PTPRC",
  "CD3D",
  "NCAM1")

p <- plotEmbedding(
  ArchRProj = subset, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  imputeWeights = NULL
)

p.genes <- lapply(p, function(x){
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
do.call(cowplot::plot_grid, c(list(ncol = 3),p.genes))


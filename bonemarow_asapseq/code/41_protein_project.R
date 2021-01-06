library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(Matrix)

# Load the reference dataset
#InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
DefaultAssay(bm) <- 'ADT'

# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)
bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model = TRUE)
ref.obj <- bm; rm(bm)

# Import the HQ cells
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Run Yuhan's transfer code
adt.matrix <- data.matrix(readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode])

rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD56(NCAM)")] <- "CD56"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD4-2")] <- "CD4"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD3-1")] <- "CD3"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD57_Recombinant" )] <- "CD57"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD38-1" )] <- "CD38"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "HLA-DR" )] <- "HLA.DR"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD278")] <-  "CD278-ICOS"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD11a/CD18")] <- "CD11a"
rownames(adt.matrix)[ which(  rownames(adt.matrix) ==  "CD127")] <- "CD127-IL7Ra"
adt.feature <- intersect(  rownames(ref.obj[["ADT"]]),  rownames(adt.matrix) )
setdiff(  rownames(ref.obj[["ADT"]]) ,  adt.feature)

# get supervised adt pca from reference
DefaultAssay(ref.obj) <- "ADT"
ref.obj <- RunSPCA(object = ref.obj, features = adt.feature, assay = "ADT", reduction.name = "sapca", reduction.key = "SAPC_", graph = "wsnn")
ref.obj <- ref.obj[adt.feature,]
ref.obj <- FindNeighbors(ref.obj)

# create another assay with matched ADTs
query <- CreateSeuratObject(counts =  adt.matrix[adt.feature,], assay = "ADT2")
DefaultAssay(query) <- "ADT2"
query <- NormalizeData(query, normalization.method = "CLR", margin = 2)
VariableFeatures(query ) <- adt.feature
query <- ScaleData(query) %>% RunPCA() %>% FindNeighbors()
anchors.transfer <- FindTransferAnchors(
  reference = ref.obj,
  query = query,
  reference.assay = "ADT",
  query.assay = 'ADT2', 
  reference.reduction = 'sapca',
  features = adt.feature,
  dims = 1:18,
  nn.method = "annoy",
  k.filter = NA,
  verbose = TRUE
)

### Find RPCA anchors, and replace the anchors in anchors.transfer
query.inte <- query
query.inte <- RunPCA(query.inte, features = adt.feature)
ref.obj.inte <- ref.obj
ref.obj.inte[["pca"]] <- ref.obj.inte[["sapca"]]
anchor.inte <- FindIntegrationAnchors(object.list = list(ref.obj.inte, query.inte),
                                      reference = c(1),
                                      anchor.features = adt.feature,
                                      dims = 1:18, 
                                      reduction = "rpca", 
                                      k.filter = NA)
# add integration anchors to transfer anchors
anchors.transfer@anchors <- anchor.inte@anchors[anchor.inte@anchors$dataset1 == 1 , 1:3]

query <- MapQuery(
  reference = ref.obj,
  query = query ,
  anchorset = anchors.transfer,
  refdata = list(celltype.l1 = "celltype.l1", 
                 celltype.l2 ="celltype.l2"),
  reference.reduction = "sapca",
  reduction.model = "wnn.umap"
)

table(query@meta.data$predicted.celltype.l1)
table(query@meta.data$predicted.celltype.l2)

predictions_l1 <- TransferData(anchorset = anchors.transfer, refdata = ref.obj$celltype.l1, dims = 1:18)
table(predictions_l1$predicted.id)


# Add ArchR clusters
# Set up color annotation / reannotion convention
cluster_reanno <- c("E1", "E2", "T1", "T2", "T3", "T4", "T5",
                    "Proj", "pDC", "B1", "B2", "B3", "B4", "Baso",
                    "M1", "M2", "M3", "M4", "M5", "M6", "M7")
names(cluster_reanno) <- paste0("C", as.character(1:21))
mdf$archr_cluster <- cluster_reanno[as.character(mdf$Clusters)]
mdf$predictions_l1 <- query@meta.data$predicted.celltype.l1
mdf$predictions_l2 <- query@meta.data$predicted.celltype.l2
mdf$confidence_l1 <-  query@meta.data$predicted.celltype.l1.score
mdf$confidence_l2 <-  query@meta.data$predicted.celltype.l2.score

ggplot(mdf, aes(x = archr_cluster, y = confidence_l1)) +
  geom_boxplot()

saveRDS(mdf, file = "../output/Seurat_Proteinprojections_allMeta.rds")

color_vec <- c("#FB6A4A", "#A50F15",  jdb_palette("brewer_marine")[c(3,5,6,7,8)],
               "green4", "#094078", "#581F99", "#6D49AB", "#8880C4", "#9DA8D7", "#CD8500",
               jdb_palette("brewer_fire")[c(4,2,5,6,7,8,9)])
names(color_vec) <- unname(cluster_reanno)

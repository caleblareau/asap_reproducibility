library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)

gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

scqc <- fread("../data/three_comparison/broad_pbmcs_aggr_singlecell.csv.gz") %>% dplyr::filter(cell_id != "None")
cellsA <- scqc[[1]][grepl("1", scqc[[1]])]
cellsB <- scqc[[1]][grepl("2", scqc[[1]])]
cellsC <- scqc[[1]][grepl("3", scqc[[1]])]

# create a gene by cell matrix
gene.activitiesA <- FeatureMatrix(
  fragments = "../../../asap_large_data_files/broad_experiment_pbmcs/cellranger_individual/ASAP_A_v12_hg38_fragments.tsv.gz",
  features = genebodyandpromoter.coords,
  cells = cellsA,
  chunk = 20
)
gene.activitiesB <- FeatureMatrix(
  fragments = "../../../asap_large_data_files/broad_experiment_pbmcs/cellranger_individual/ASAP_B_v12_hg38_fragments.tsv.gz",
  features = genebodyandpromoter.coords,
  cells = cellsB,
  chunk = 20
)
gene.activitiesC <- FeatureMatrix(
  fragments = "../../../asap_large_data_files/broad_experiment_pbmcs/cellranger_individual/ASAP_C_v12_hg38_fragments.tsv.gz",
  features = genebodyandpromoter.coords,
  cells = cellsC,
  chunk = 20
)
gene.activies <- cbind(gene.activitiesA, gene.activitiesB, gene.activitiesC)
# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_name
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

saveRDS(gene.activities, file = "../../../asap_large_data_files/broad_experiment_pbmcs/output/ASAP_BROAD_aggr_29April2020.rds")

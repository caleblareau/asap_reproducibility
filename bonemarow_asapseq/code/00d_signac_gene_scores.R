library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)

# Get coordinates
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

# Get fragment files
frag_file <- "../../../asap_large_data_files/bonemarrow_data/input/asap_marrow_hg38_fragments.tsv.gz"
bcs <- fread("../data/barcodes/step3_ADThq.tsv")[[1]]

process_ga <- function(frag_file, cells){
  # create a gene by cell matrix
  gene.activities <- FeatureMatrix(
    fragments = frag_file,
    features = genebodyandpromoter.coords,
    cells = cells,
    chunk = 20
  )
  
  gene.key <- genebodyandpromoter.coords$gene_name
  names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
  rownames(gene.activities) <- gene.key[rownames(gene.activities)]
  gene.activities
}

bmga <- process_ga(frag_file, bcs)

saveRDS(bmga, file = "../../../asap_large_data_files/bonemarrow_data/output/signac_marrow_gene_scores.rds")

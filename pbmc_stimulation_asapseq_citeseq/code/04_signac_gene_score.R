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
frag_files <- paste0("../../../asap_large_data_files/pbmc_stim_data/input/",c("control", "stim"),"_fragments.tsv.gz")
cdf <- readRDS("../output/adt_mat/ASAP_embedding_CLRadt.rds")

# Identify cells to score from ArchR data
control_cells <- gsub("Control#", "", rownames(cdf)[cdf$sample == "Control"])
stim_cells <- gsub("Stim#", "", rownames(cdf)[cdf$sample == "Stim"])

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

ga_control <- process_ga(frag_files[1], control_cells)
ga_stim <- process_ga(frag_files[2], stim_cells)

saveRDS(ga_control, file = "../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_control_ga.rds")
saveRDS(ga_stim, file = "../../../asap_large_data_files/pbmc_stim_data/output/signac_genes_scores/signac_stim_ga.rds")

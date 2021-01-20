library(Seurat)
library(data.table)
library(Signac)

# Import protein
import_kite_counts <- function(lib_id, barcode_idx = ""){
  mtx <- fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]], "-", barcode_idx)
  colnames(matx) <- paste0(fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

lll_protein <- cbind(import_kite_counts("LLL_ctrl", "1"),import_kite_counts("LLL_stim", "2"))

# Import RNAseq
# Import scRNA-seq
process_scRNA  <- function(raw){
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 30 & n_total_rna > 1000 & n_feature_rna > 500]
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), qc_cells]
  raw
}

lll_stim_rna <- Read10X_h5("../data/LLL_STIM_feature_bc_matrix.h5")$`Gene Expression`
colnames(lll_stim_rna) <- gsub("-1", "-2", colnames(lll_stim_rna))
lll_stim_rna <- process_scRNA(lll_stim_rna)

lll_ctrl_rna <- process_scRNA(Read10X_h5("../data/LLL_CTRL_feature_bc_matrix.h5")$`Gene Expression`)
pbmc_lll <- CreateSeuratObject(counts = cbind(lll_ctrl_rna,lll_stim_rna))

# Import ATAC-seq
peaks <- Read10X_h5("../../../asap_large_data_files/multiome_pbmc_stim/input/filtered_peak_bc_matrix.h5")
lll_peaks <- peaks[,substr(colnames(peaks), 18,18) %in% c("1", "2")][,colnames(pbmc_lll)]

chrom_assay <- CreateChromatinAssay(
  counts = lll_peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../../asap_large_data_files/multiome_pbmc_stim/input/fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)

pbmc_lll[["peaks"]] <- chrom_assay
pbmc_lll[["ADT"]] <- CreateAssayObject(lll_protein[,colnames(pbmc_lll)])

# RNA Seurat stuff
pbmc_lll <- FindVariableFeatures(pbmc_lll)
pbmc_lll <- NormalizeData(pbmc_lll)
pbmc_lll <- ScaleData(pbmc_lll)

# Do it for ADT
DefaultAssay(pbmc_lll) <- "ADT"
pbmc_lll <- FindVariableFeatures(pbmc_lll)
pbmc_lll <- NormalizeData(pbmc_lll)
pbmc_lll <- ScaleData(pbmc_lll)

saveRDS(pbmc_lll, "../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")


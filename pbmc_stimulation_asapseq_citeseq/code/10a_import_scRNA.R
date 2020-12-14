# Import scRNA-seq
import_scRNAseq_cite_stim <- function(dir_base){
  
  data.dir <- paste0("../../pbmc_stimulation_citeseq/data/rnaseq/", dir_base)
  raw <- Read10X(data.dir = data.dir)
  colnames(raw) <- paste0(substr(colnames(raw), 1, 16), "-1")
  
  # import scrublet results
  singlets <- fread(paste0("../../pbmc_stimulation_citeseq/data/rnaseq/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(score < 0.25) %>% pull(barcode) # the original called threshold seemed too conservative; this is a better estimate for these libraries
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 10 & n_total_rna > 1000 & n_feature_rna > 500]
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), intersect(singlets, qc_cells)]
  raw <- CreateSeuratObject(counts = raw, project = "RNA")
  raw <- FindVariableFeatures(raw)
  raw <- NormalizeData(raw)
  raw <- ScaleData(raw)
  raw
}
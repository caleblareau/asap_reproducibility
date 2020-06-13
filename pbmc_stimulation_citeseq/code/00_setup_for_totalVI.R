library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)

# Import scRNA-seq
import_scRNAseq <- function(dir_base){
  
  data.dir <- paste0("../data/rnaseq/", dir_base)
  raw <- Read10X(data.dir = data.dir)
  
  # import scrublet results
  singlets <- fread(paste0("../data/rnaseq/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(score < 0.2) %>% pull(barcode) # the original called threshold seemed too conservative; this is a better estimate for these libraries
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), singlets]
  raw
}

# Import ADT
import_kite_counts <- function(library){
  mtx <- fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

export_for_totalVI <- function(adt_list, rna_list, out_dir){
  
  stopifnot(length(adt_list) == length(rna_list))
  
  # Remove isotype controls from ADT
  protein_feature_names <- rownames(adt_list[[1]])[!grepl("sotypeCtrl", rownames(adt_list[[1]]))]
  
  # Create one giant matrix of protein count values
  adt_mat <- do.call("cbind", lapply(1:length(adt_list), function(idx){
    m <- adt_list[[idx]][protein_feature_names,]
    colnames(m) <- paste0(substr(colnames(m), 1, 16), "-", as.character(idx))
    m
  }))
  
  # Create large binary ATAC matrix
  rna_mat <- do.call("cbind", lapply(1:length(rna_list), function(idx){
    m <- rna_list[[idx]]
    colnames(m) <- paste0(substr(colnames(m), 1, 16), "-", as.character(idx))
    m
  }))
  
  # Make sure they have the same cell names
  stopifnot(colnames(rna_mat) == colnames(adt_mat))

  # Create feature names
  rna_feature_names <- rownames(rna_list[[1]])
  
  
  # Prepare for export
  features_df <- data.frame(
    x1 = c(rna_feature_names, protein_feature_names),
    x2 = c(rna_feature_names, protein_feature_names),
    x3 = c(rep("Gene Expression", length(rna_feature_names)),
           rep("Antibody Capture", length(protein_feature_names)))
  )
  
  out_mat <- rbind(rna_mat, adt_mat)
  
  out_barcodes <- data.frame(barcode = colnames(out_mat))
  
  # Make relevant landing place for the data
  dir.create(out_dir)
  dir.create(paste0(out_dir, "/filtered_feature_bc_matrix"))
  
  # Export matrix
  out_mat_file = paste0(out_dir, "/filtered_feature_bc_matrix/matrix.mtx")
  writeMM(out_mat, out_mat_file)
  system(paste0("gzip ", out_mat_file))
  
  # Export barcodes
  out_barcode_file <- paste0(out_dir, "/filtered_feature_bc_matrix/barcodes.tsv")
  write.table(out_barcodes, file = out_barcode_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  system(paste0("gzip ", out_barcode_file))
  
  # Export features
  out_features_file <- paste0(out_dir, "/filtered_feature_bc_matrix/features.tsv")
  write.table(features_df, file = out_features_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
  system(paste0("gzip ", out_features_file))
  
  # Nomial return value
  out_dir
}
# Import matrices
ctrl_scRNA <- import_scRNAseq("ctrl")
stim_scRNA <- import_scRNAseq("stim")
ctrl_ADT <- import_kite_counts("ctrl")[,gsub("-1", "", colnames(ctrl_scRNA))]
stim_ADT <- import_kite_counts("stim")[,gsub("-1", "", colnames(stim_scRNA))]

booctrl <- (colSums(ctrl_ADT) > 25000)
boostim <- (colSums(stim_ADT) > 25000)

# Aggregate
export_for_totalVI(list(ctrl_ADT[,booctrl], stim_ADT[,boostim]), list(ctrl_scRNA[,booctrl], stim_scRNA[,boostim]),
                   "../output/pbmcstim_citeseq_combined")




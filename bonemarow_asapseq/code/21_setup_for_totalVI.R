library(SummarizedExperiment)
library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggrastr)
library(BuenColors)

ga_ss <- readRDS("../../../asap_large_data_files/bonemarrow_data/output/signac_marrow_gene_scores.rds")
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,colnames(ga_ss)]
rownames(adt_ss) <- gsub("_", "-", rownames(adt_ss))
adt_ss <- adt_ss[!grepl("sotype", rownames(adt_ss)),]

export_for_scVI <- function(adt_list, ga_list, out_dir, ngenes = 3000){
  
  stopifnot(length(adt_list) == length(ga_list))
  
  # Create one giant matrix of protein count values
  adt_mat <- do.call("cbind", lapply(1:length(adt_list), function(idx){
    m <- adt_list[[idx]]
    colnames(m) <- paste0(substr(colnames(m), 1, 16), "-", as.character(idx))
    m
  }))
  
  # Create large binary ATAC matrix
  atac_mat <- do.call("cbind", lapply(1:length(ga_list), function(idx){
    m <- ga_list[[idx]]
    colnames(m) <- paste0(substr(colnames(m), 1, 16), "-", as.character(idx))
    m
  }))
  
  # Make sure they have the same cell names
  stopifnot(colnames(atac_mat) == colnames(adt_mat))
  
  # Filter for top ngenes
  rs <- rowSums(atac_mat)
  thresh <- sort(rs, decreasing = TRUE)[ngenes]
  atac_mat <- atac_mat[rs >= thresh,]
  
  # Create feature names
  all_feature_names <- rownames(ga_list[[1]])
  atac_feature_names <- all_feature_names[rs >= thresh]
  protein_feature_names <- rownames(adt_list[[1]])
  
  # Prepare for export
  features_df <- data.frame(
    x1 = c(atac_feature_names, protein_feature_names),
    x2 = c(atac_feature_names, protein_feature_names),
    x3 = c(rep("Gene Expression", length(atac_feature_names)),
           rep("Antibody Capture", length(protein_feature_names)))
  )
  
  out_mat <- rbind(atac_mat, adt_mat)
  
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

export_for_scVI(list(adt_ss), list(ga_ss), "../output/for_total_vi_asap_marrow")


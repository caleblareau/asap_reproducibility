library(data.table)
library(Matrix)

import_kite_counts_NYGC <- function(path, library, gz = FALSE, string_col = "", string_row = ""){
  
  xx <- ifelse(gz, ".gz", "")
  mtx <- fread(paste0(path, "featurecounts",library,".mtx",xx), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".barcodes.txt",xx), header = FALSE)[[1]], "-1")
  colnames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".genes.txt",xx), header = FALSE)[[1]])
  rownames(matx) <- paste0(string_row, rownames(matx))
  colnames(matx) <- paste0(string_col, colnames(matx))
  
  return(matx)
}
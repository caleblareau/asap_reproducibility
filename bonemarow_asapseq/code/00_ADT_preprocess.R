library(data.table)
library(Matrix)

# Do a small first pass filter to ID
# The super GEMs and the very low counts

# Import ADT
import_kite_counts_bm <- function(){
  mtx <- fread(paste0("../data/adt/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../data/adt/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../data/adt/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

bm <- import_kite_counts_bm()
bm <- bm[,colSums(bm) >= 100]

# Filter out very high / isotype control guys
bm_meta <- data.frame(
  total = colSums(bm),
  control_tags = colSums(bm[grepl("sotypeCtrl", rownames(bm)),]) + 1
)

ggplot(bm_meta, aes(total, control_tags)) +
  geom_point() + scale_x_log10() + scale_y_log10()

# Get rid of the super expressers
remove_bm <- bm_meta$total > 10000 | bm_meta$control_tags > 50
sum(remove_bm)
bm <- bm[,!remove_bm]

saveRDS(bm, file = "../output/adt_mat/ASAP_bonemarrow.rds")



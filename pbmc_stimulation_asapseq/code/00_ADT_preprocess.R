library(data.table)
library(Matrix)

# Do a small first pass filter to ID
# The super GEMs and the very low counts

# Import ADT
import_kite_counts <- function(library){
  mtx <- fread(paste0("../data/adt/ASAP_",library,"_ADT/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../data/adt/ASAP_",library,"_ADT/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../data/adt/ASAP_",library,"_ADT/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

ctrl <- import_kite_counts("ctrl")
stim <- import_kite_counts("stim")

# Filter very low counts
ctrl <- ctrl[,colSums(ctrl) >= 100]
stim <- stim[,colSums(stim) >= 100]

# Filter out very high / isotype control guys
ctrl_meta <- data.frame(
  total = colSums(ctrl),
  control_tags = colSums(ctrl[grepl("sotypeCtrl", rownames(ctrl)),]) + 1
)

stim_meta <- data.frame(
  total = colSums(stim),
  control_tags = colSums(stim[grepl("sotypeCtrl", rownames(stim)),]) + 1
)

ggplot(ctrl_meta, aes(total, control_tags)) +
  geom_point() + scale_x_log10() + scale_y_log10()

remove_stim <- stim_meta$total > 25000 | stim_meta$control_tags > 75
remove_ctrl <- ctrl_meta$total > 25000 | ctrl_meta$control_tags > 75

sum(remove_stim)
sum(remove_ctrl)

stim <- stim[,!remove_stim]
ctrl <- ctrl[,!remove_ctrl]

saveRDS(stim, file = "../output/adt_mat/ASAP_stim.rds")
saveRDS(ctrl, file = "../output/adt_mat/ASAP_ctrl.rds")

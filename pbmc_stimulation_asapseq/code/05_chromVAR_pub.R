library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(chromVAR)

# function to get counts
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)
  
  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found
  
  # Make a polished colData
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}

# Import PBMC to a summarized experiment
bc_control <- as.character(read.table("../data/cellranger/control_barcodes.tsv")[,1])
bc_stim <- as.character(read.table("../data/cellranger/stim_barcodes.tsv")[,1])
peaks <- diffloop::bedToGRanges("../data/published/CD4_stimulation.bed")

control_SE <- getCountsFromFrags("../../../asap_large_data_files/pbmc_stim_data/input/control_fragments.tsv.gz", peaks, bc_control)
stim_SE <- getCountsFromFrags("../../../asap_large_data_files/pbmc_stim_data/input/stim_fragments.tsv.gz", peaks, bc_stim)

files <- c("../data/published/CD3_CD28_gained_peaks.bed"); names <- c("gained")
mat <- sapply(1:length(files), function(i){
  v <- as.numeric(1:length(peaks) %in% queryHits(findOverlaps(peaks, diffloop::bedToGRanges(files[i]))))
  return(v)
})
colnames(mat) <- names
mmat <- Matrix(mat)
remove(mat)

library(BSgenome.Hsapiens.UCSC.hg38)
control_SE <- addGCBias(control_SE, BSgenome.Hsapiens.UCSC.hg38)
stim_SE <- addGCBias(stim_SE, BSgenome.Hsapiens.UCSC.hg38)
boo_peak_stim <- rowSums(assays(stim_SE)[["counts"]]) >0
boo_peak <- rowSums(assays(control_SE)[["counts"]]) >0
stim_SE <- stim_SE[boo_peak &boo_peak_stim, ]
control_SE <- control_SE[boo_peak &boo_peak_stim, ]
bg <- getBackgroundPeaks(control_SE, niterations = 100)

control_score <- data.matrix(assays(computeDeviations(control_SE, mmat[boo_peak & boo_peak_stim, ,drop = FALSE], background_peaks = bg))[["z"]])
stim_score <- data.matrix(assays(computeDeviations(stim_SE[boo_peak_stim], mmat[boo_peak &boo_peak_stim,, drop = FALSE], background_peaks = bg))[["z"]])
qplot(stim_score[1,])
qplot(control_score[1,])


library(dplyr)
library(data.table)
library(stringr)


# Assign the bus set to the corresponding 
process_tag_reference <- function(ref_file){
 
  ref_df <- fread(ref_file, header = FALSE)
  ref_df2 <- data.frame(
    seq = ref_df[[1]][c(FALSE,TRUE)],
    name = gsub(">", "", ref_df[[1]][c(TRUE,FALSE)])
  )
  ref_df2$idx <- as.character(1:dim(ref_df2)[1])
  ref_df2$name_short <- str_split_fixed(ref_df2$name, "-", 2)[,1]
  
  # Some of these Abs have an extra -# that is necessary
  doubles <- data.frame(table(ref_df2$name_short)) %>% filter(Freq > 50) %>% pull(Var1) %>% as.character()
  ref_df2$extra_idx <- str_split_fixed(ref_df2$name, "-", 3)[,2]
  ref_df2$use <- ifelse(ref_df2$name_short %in% doubles, paste0(ref_df2$name_short,"_", ref_df2$extra_idx),ref_df2$name_short)
  
  vec <- as.character(ref_df2$use); names(vec) <- as.character(ref_df2$idx)
  vec
}

# 
pbmc_vec <- process_tag_reference("references/PBMC_FeaturesMismatch.fa")
BM_vec <- process_tag_reference("references/BM_FeaturesMismatch.fa")

# Function translated from java version: https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java
# Not vectorized!!! 
estimateLibrarySize <- function(nTotal, nUnique){
  
  f <- function(x, c, n) {
    return(c / x - 1 + exp(-n / x))
  }
  
  m = 1
  M = 100
  
  nDuplicates <- (nTotal - nUnique) + 1 # charity to handle only unique reads observed
  
  if (nUnique > nTotal | (f(m * nUnique, nUnique, nTotal) < 0) | nUnique < 0 | nTotal < 0 | nDuplicates < 0) {
    message("Library size returns 0 -- invalid inputs; check this cell more closely")
    return(0)
  }
  
  while (f(M * nUnique, nUnique, nTotal) > 0) {
    M <- M*10.0
  }
  
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nUnique, nUnique, nTotal);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }
  
  return(round(nUnique * (m + M) / 2.0))
}

# Function to process the per-barcode summary statistics
estimate_tag_rates <- function(bus_file, barcode_file, vec){
  bcs <- substr(fread(barcode_file, header = FALSE)[[1]], 1, 16)
  dff <- fread(bus_file, col.names = c("barcode", "umi", "set", "idx")) %>%
    mutate(Ab = vec[as.character(set)], umi_ab = paste0(Ab,"_",umi)) %>%
    filter(barcode %in% bcs) %>% 
    group_by(barcode) %>%
    summarise(count = n_distinct(umi_ab), n_total = n()) %>%
    mutate(n_duplicated = n_total - count) 
  dff$complexity <- sapply(1:(dim(dff)[1]), function(i){
    estimateLibrarySize(dff$n_total[i],dff$count[i])
  })
  dff
}

ASAP_ctrl <- estimate_tag_rates("ASAP_ctrl_featurecounts.tsv", "../../CONTROL_v12_hg38-mtMask_FC14k/outs/filtered_tf_bc_matrix/barcodes.tsv", pbmc_vec)
ASAP_stim <- estimate_tag_rates("ASAP_stim_featurecounts.tsv", "../../STIM_v12_hg38-mtMask_FC14k/outs/filtered_tf_bc_matrix/barcodes.tsv", pbmc_vec)
CITE_ctrl <- estimate_tag_rates("CITE_ctrl_featurecounts.tsv", "../../../../asap_paper/asap_reproducibility/pbmc_stimulation_citeseq/data/rnaseq/ctrl/barcodes.tsv.gz", pbmc_vec)
CITE_stim <- estimate_tag_rates("CITE_stim_featurecounts.tsv", "../../../../asap_paper/asap_reproducibility/pbmc_stimulation_citeseq/data/rnaseq/stim/barcodes.tsv.gz", pbmc_vec)
ASAP_marrow <- estimate_tag_rates("BM_featurecounts.tsv", "../../../BM_asap_v12_hg38-mtMask_mgatk/final/BM_asap_v12_hg38-mtMask_mgatk.depthTable.txt", BM_vec)

summary(ASAP_ctrl[ASAP_ctrl$complexity < 100000, "count"])
summary(CITE_ctrl[CITE_ctrl$complexity < 100000, "count"])
summary(ASAP_stim[ASAP_stim$complexity < 100000, "count"])
summary(CITE_stim[CITE_stim$complexity < 100000, "count"])
summary(ASAP_marrow[ASAP_marrow$complexity < 100000, "count"])

compute_duplicate_percentage <- function(df){
  df %>% filter(complexity < 100000) %>%
    mutate(percent_duplicates = n_duplicated/n_total*100) %>%
    pull(percent_duplicates) %>% summary
}
compute_duplicate_percentage(ASAP_ctrl)
compute_duplicate_percentage(CITE_ctrl)
compute_duplicate_percentage(ASAP_ctrl)
compute_duplicate_percentage(ASAP_stim)
compute_duplicate_percentage(ASAP_marrow)




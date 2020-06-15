library(data.table)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)

"%ni%" <- Negate("%in%")

cdf <- readRDS("../output/adt_mat/ASAP_embedding_CLRadt.rds")
control_frags <- fread("../../../asap_large_data_files/pbmc_stim_data/input/control_fragments.tsv.gz", header = FALSE) %>% data.frame()
stim_frags <- fread("../../../asap_large_data_files/pbmc_stim_data/input/stim_fragments.tsv.gz", header = FALSE) %>% data.frame()

create_bigwig <- function(frags, clusters, name, condition){
  
  cdf$cell_id <- rownames(cdf)
  possible_ids <- cdf %>% filter(cluster %in% clusters) %>%
    filter(sample == condition)  %>% 
    pull(cell_id) %>% as.character()
  cr_ids <- stringr::str_split_fixed(possible_ids, "#", 2)[,2]
  cluster_gr <- frags %>% filter(V4 %in% cr_ids) %>%
    setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
    makeGRangesFromDataFrame()
  
  reads_coverage <- coverage(cluster_gr)/length(cluster_gr)*1000000
  export.bw(reads_coverage, con = paste0("../../../asap_large_data_files/pbmc_stim_data/output/pseudobulk_bigwigs/",name,"_",condition,".bw"))
  name
}

# Create bigwigs for major celltypes
create_bigwig(control_frags, c("C5"), "Bcell", "Control")
create_bigwig(stim_frags   , c("C5"), "Bcell", "Stim")

create_bigwig(control_frags, c("C8", "C9"), "NK", "Control")
create_bigwig(stim_frags   , c("C8", "C9"), "NK", "Stim")

create_bigwig(control_frags, c("C7", "C11"), "CD8T", "Control")
create_bigwig(stim_frags   , c("C7", "C11"), "CD8T", "Stim")

create_bigwig(control_frags, c("C10", "C12", "C13", "C14", "C15"), "CD4T", "Control")
create_bigwig(stim_frags   , c("C10", "C12", "C13", "C14", "C15"), "CD4T", "Stim")



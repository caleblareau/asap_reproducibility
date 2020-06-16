library(data.table)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)

"%ni%" <- Negate("%in%")

cdf <- readRDS("../output/adt_mat/ASAP_embedding_CLRadt.rds")
control_frags <- fread("../../../asap_large_data_files/pbmc_stim_data/input/control_fragments.tsv.gz", header = FALSE) %>% data.frame()
stim_frags <- fread("../../../asap_large_data_files/pbmc_stim_data/input/stim_fragments.tsv.gz", header = FALSE) %>% data.frame()

frags <- rbind(control_frags, stim_frags)
cdf$cell_id <- rownames(cdf)

# Bright/dim analysis
clusters <- c("C8", "C9")
possible_ids_bright <- cdf %>% filter(cluster %in% clusters) %>%
  filter(CD56.NCAM. >= 3) %>% 
  pull(cell_id) %>% as.character()

possible_ids_dim <- cdf %>% filter(cluster %in% clusters) %>%
  filter(CD56.NCAM. < 3) %>% 
  pull(cell_id) %>% as.character() 

# Modify archr names
dim_ids <- stringr::str_split_fixed(possible_ids_dim, "#", 2)[,2]
bright_ids <- stringr::str_split_fixed(possible_ids_bright, "#", 2)[,2]

# Make GRanges
cluster_gr_dim <- frags %>% filter(V4 %in% dim_ids) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

cluster_gr_bright <- frags %>% filter(V4 %in% bright_ids) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn")) %>%
  makeGRangesFromDataFrame()

reads_coverage_bright <- coverage(cluster_gr_bright)/length(cluster_gr_bright)*1000000
export.bw(reads_coverage_bright, con = paste0("../../../asap_large_data_files/pbmc_stim_data/output/pseudobulk_bigwigs/NK_CD56_bright.bw"))

reads_coverage_dim <- coverage(cluster_gr_dim)/length(cluster_gr_dim)*1000000
export.bw(reads_coverage_dim, con = paste0("../../../asap_large_data_files/pbmc_stim_data/output/pseudobulk_bigwigs/NK_CD56_dim.bw"))

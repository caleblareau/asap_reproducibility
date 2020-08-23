library(data.table)
library(BuenColors)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

load("../output/AvsB_LLL-Omni.rda")
plot_df <-  data.frame(archr_df, adt_mat_CLR)

ggplot(plot_df, aes(x =IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = TSA_CD56) ) +
  geom_point() +
  scale_color_gradientn(colors = jdb_palette("solar_extra"))

ggplot(plot_df, aes(x =IterativeLSI.UMAP_Dimension_1, y = IterativeLSI.UMAP_Dimension_2, color = Clusters) ) +
  geom_point() + scale_color_manual(values = jdb_palette("corona"))

plot_df %>% filter(Clusters == "C3") %>%
  ggplot(aes(x = TSA_CD56)) + 
  geom_histogram(bins = 50)

bright_barcodes <- plot_df %>% filter(Clusters == "C3") %>%
  filter(TSA_CD56 > 0.75) %>% pull(barcode)

dim_barcodes <- plot_df %>% filter(Clusters == "C3") %>%
  filter(TSA_CD56 < 0.5)%>% pull(barcode)

LLL_frags <- fread("../../../asap_large_data_files/nygc_pbmc/input/AvsB_ATAC_LLL_fragments.tsv.gz")
OMNI_frags <- fread("../../../asap_large_data_files/nygc_pbmc/input/AvsB_ATAC_OMNI_fragments.tsv.gz")

LLL_frags$barcode <- paste0("AvsB_LLL#", LLL_frags$V4)
OMNI_frags$barcode <- paste0("AvsB_OMNI#", OMNI_frags$V4)
frags <- rbind(LLL_frags, OMNI_frags)

# Make GRanges
cluster_gr_dim <- frags %>% filter(barcode %in% dim_barcodes) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn", "barcode")) %>%
  makeGRangesFromDataFrame()

cluster_gr_bright <- frags %>% filter(barcode %in% bright_barcodes) %>%
  setnames(c("chr", "start", "end", "V4", "PCRn", "barcode")) %>%
  makeGRangesFromDataFrame()

reads_coverage_bright <- coverage(cluster_gr_bright)/length(cluster_gr_bright)*1000000
export.bw(reads_coverage_bright, con = paste0("../output/NK_CD56_bright.bw"))

reads_coverage_dim <- coverage(cluster_gr_dim)/length(cluster_gr_dim)*1000000
export.bw(reads_coverage_dim, con = paste0("../output/NK_CD56_dim.bw"))

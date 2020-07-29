library(data.table)

mdf <- data.frame(readRDS("../output/ArchR_main_metadata.rds"))
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))
basophil_barcodes <- mdf %>% filter(Clusters == "C14") %>% pull(barcode)
frags <- fread("../../../asap_large_data_files/bonemarrow_data/input/asap_marrow_hg19_fragments.tsv.gz")
frags_b <- frags %>% filter(V4 %in% basophil_barcodes)
write.table(frags_b[,c(1,2,3,4)], file = "../output/basophils/hg19_basophil_frags.bed",
            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

frags_b_gr <- makeGRangesFromDataFrame(frags_b, 
                                       seqnames.field = "V1", start.field = "V2", end.field = "V3")
cov <- coverage(frags_b_gr)/length(frags_b_gr)*1000000
rtracklayer::export.bw(cov, con = paste0("../output/basophils/hg19_basophil_asap_marrow.bw"))

# Call peaks

# Now import gwas / prev peaks
baso_gwas<- fread("../data/pub/baso_count.bext") %>% data.frame() %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "V1", start.field = "V2", end.field = "V3")
prev_gr <- fread("../data/pub/29August2017_EJCsamples_allReads_500bp.bed", col.names = c("chr", "start", "end")) %>% data.frame() %>% makeGRangesFromDataFrame()
baso_peaks <- fread("../output/basophils/NA_peaks.narrowPeak") %>% data.frame() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqnames.field = "V1", start.field = "V2", end.field = "V3")

ov_prev <- findOverlaps(baso_gwas, prev_gr)
ov_new <- findOverlaps(baso_gwas, baso_peaks)

queryHits(ov_new) %in% queryHits(ov_prev)
data.frame(baso_gwas)[1:length(baso_gwas) %in% queryHits(ov_new),] %>% 
  arrange(desc(V5))

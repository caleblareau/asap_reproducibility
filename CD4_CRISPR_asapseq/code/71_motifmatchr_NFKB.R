options(stringsAsFactors=F)

#******************************************
### args        = commandArgs(trailingOnly = T)
### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP="Perturb_CD4_stim"

#******************************************
suppressMessages(library(motifmatchr))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(GenomicRanges))
suppressMessages(library(magrittr))
suppressMessages(library(chromVARmotifs))
suppressMessages(library(tidyverse))

dir.create("../output/motif_NFKB2",   showWarnings = FALSE, recursive = TRUE)

#******************************************
data("human_pwms_v1")
# Pull out the motifs that you want
NFKB2 = human_pwms_v1[names(human_pwms_v1)[grepl("NFKB2", names(human_pwms_v1))]][[1]]

#Import peaks
peaks_IL2RA  = read.table("../output/each_sgRNA/MACS2_peak_call_merge.bed", header = FALSE, col.names = c("chr", "start", "end")) %>%
                filter(chr=="chr10",start>6000000,end<6100000) %>% 
                makeGRangesFromDataFrame(.)
motif_res     = matchMotifs(NFKB2, peaks_IL2RA, BSgenome.Hsapiens.UCSC.hg38, out = "positions")
#### GRangesList object of length 1:
#### [[1]]
#### GRanges object with 8 ranges and 1 metadata column:
####       seqnames          ranges strand |            score
####          <Rle>       <IRanges>  <Rle> |        <numeric>
####   [1]    chr10 6045258-6045270      + | 7.01299700991556
####   [2]    chr10 6045258-6045270      - | 7.30396707098777
####   [3]    chr10 6062567-6062579      + |  7.9547807077044
####   [4]    chr10 6062568-6062580      + | 8.65783924621558
####   [5]    chr10 6062567-6062579      - |  8.1241572267519
####   [6]    chr10 6062568-6062580      - | 8.97386070041671
####   [7]    chr10 6088754-6088766      + | 9.75890208505077
####   [8]    chr10 6088754-6088766      - | 9.20092660589381
####   -------
####   seqinfo: 1 sequence from an unspecified genome; no seqlengths

motif_res_bed = motif_res %>% as.data.frame(.) %>% select(seqnames,start,end) %>% unique()
write.table(motif_res_bed,"output/motif_NFKB2/motif_NFKB2.bed",row.names=F, col.names=F, sep="\t", append=F, quote=F)

# Look at first motif
x <- data.matrix(NFKB2@profileMatrix)
m <- round(-0.258 * (0.0310078- exp(x)), 4)

pdf("../output/motif_NFKB2/motif_NFKB2_seqLogo.pdf",h=5,w=5)
 seqLogo(m)
dev.off()


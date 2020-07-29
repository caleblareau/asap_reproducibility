library(data.table)
library(BuenColors)

cd4_bed <- fread("CD4_stimulation.bed")
cd4_counts <- fread("CD4_stimulation.counts.tsv")

library(DESeq2)
ddsMat <- DESeqDataSetFromMatrix(countData = cd4_counts,
                                 colData = data.frame(
                                   time = c("T0", "T0", "T0", "T24", "T24", "T24")
                                 ),
                                 design = ~ time)
ddsMat <- DESeq(ddsMat)
rdf <- results(ddsMat)


#write.table(cd4_bed[which(rdf$log2FoldChange > 1 & rdf$padj < 0.01),], 
#            file = "CD3_CD28_gained_peaks.bed", 
#            row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

ggplot(data.frame(rdf), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() + pretty_plot() + L_border()

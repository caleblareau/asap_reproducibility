library(data.table)
library(dplyr)
library(Matrix)
library(dplyr)

import_cite_counts <- function(library){
  mtx <- fread(paste0("../../pbmc_stimulation_citeseq/data/adt/",library,"_featurecounts/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../../pbmc_stimulation_citeseq/data/adt/",library,"_featurecounts/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../../pbmc_stimulation_citeseq/data/adt/",library,"_featurecounts/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

import_asap_counts <- function(library){
  mtx <- fread(paste0("../../pbmc_stimulation_asapseq/data/adt/",library,"_ADT/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../../pbmc_stimulation_asapseq/data/adt/",library,"_ADT/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../../pbmc_stimulation_asapseq/data/adt/",library,"_ADT/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}


asap_control <- import_asap_counts("ASAP_ctrl")
asap_stim <- import_asap_counts("ASAP_stim")
cite_control <- import_cite_counts("ctrl")
cite_stim <- import_cite_counts("stim")

cpm <- function(mat){
  round(Matrix::rowSums(mat)/sum(mat) *1000000, 2)
}

cpm_df <- data.frame(
  tag = rownames(asap_control),
  asap_control = cpm(asap_control),
  asap_stim = cpm(asap_stim),
  cite_control = cpm(cite_control),
  cite_stim = cpm(cite_stim)
)

cpm_df$logFC_ASAP <- round(log2(cpm_df$asap_stim/cpm_df$asap_control), 3)
cpm_df$logFC_CITE <- round(log2(cpm_df$cite_stim/cpm_df$cite_control), 3)

text_plot <- ggplot(cpm_df,aes(x = logFC_ASAP, y = logFC_CITE, label = tag)) + 
  geom_text() + pretty_plot() + L_border()

dot_plot <- ggplot(cpm_df, aes(x = logFC_ASAP, y = logFC_CITE, label = tag)) + 
  geom_point(size = 0.4) + pretty_plot(fontsize = 8) + L_border() 

cowplot::ggsave2(cowplot::plot_grid(dot_plot, text_plot, nrow = 1), 
                 width = 10, height = 5.5, file = "../plots/logFC_CITE_ASAP.png")

cowplot::ggsave2(dot_plot, 
                 width = 1.5, height = 1.5, file = "../plots/logFC_CITE_ASAP_paper.pdf")

cor(cpm_df$logFC_ASAP, cpm_df$logFC_CITE)


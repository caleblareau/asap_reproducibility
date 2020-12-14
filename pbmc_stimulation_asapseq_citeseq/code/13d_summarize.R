library(data.table)
library(dplyr)

ATAC <- readRDS("../output/DGE_ATAC_Tcell_Activation.rds")
RNA <- readRDS("../output/DGE_RNA_Tcell_Activation.rds")
PROTEIN <- readRDS("../output/DGE_Protein_Tcell_Activation.rds")

ATAC_hits <- ATAC$FDR < 0.01 & (abs(ATAC$logFC) > 0.5)
RNA_hits <- RNA$FDR < 0.01 & (abs(RNA$logFC) > 0.5)
PROTEIN_hits <- PROTEIN$padj < 0.01 & (abs(PROTEIN$foldchange) > 0.5)

table(ATAC_hits)
mean(ATAC_hits)*100

table(RNA_hits)
mean(RNA_hits)*100

table(PROTEIN_hits)
mean(PROTEIN_hits)*100

length(ATAC_hits)
length(RNA_hits)
length(PROTEIN_hits)
library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
# Import features for module scoring
stim_peaks <- gsub(":", "-", readRDS("../output/top2k_stim_peaks.rds"))
stim_genes <- readRDS("../../pbmc_stimulation_asapseq_citeseq/output/DGE_RNA_Tcell_Activation.rds") %>% dplyr::filter(logFC >0) %>%
  arrange(desc(F)) %>% pull(gene) %>% head(1000)
stim_proteins <- readRDS("../../pbmc_stimulation_asapseq_citeseq/output/DGE_Protein_Tcell_Activation.rds") %>% 
  arrange(desc(foldchange)) %>% pull(marker) %>% head(20)

# Now process with module scores
pbmc_lll <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")
DefaultAssay(pbmc_lll) <- "ADT"
pbmc_lll <- AddModuleScore(pbmc_lll, features = list(stim_proteins), assay = "ADT", ctrl = 5, 
               name = "protein_MS")
DefaultAssay(pbmc_lll) <- "RNA"
pbmc_lll <- AddModuleScore(pbmc_lll, features = list(stim_genes), assay = "RNA", 
                           name = "rna_MS")
DefaultAssay(pbmc_lll) <- "peaks"
pbmc_lll <- AddChromatinModule(pbmc_lll, features = list("peaks_MS" = stim_peaks), assay = "peaks", 
                               genome = BSgenome.Hsapiens.UCSC.hg38)

pbmc_lll@meta.data$stim <- ifelse(substr(colnames(pbmc_lll), 18, 18) == 1, "Control", "Stim")

write.table(pbmc_lll@meta.data[,c("stim","protein_MS1", "rna_MS1","peaks.MS")],
            file = "../output/LLL_module_scores.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)

ggplot(pbmc_lll@meta.data, aes(x = peaks.MS, y = protein_MS1)) +
  geom_point() + facet_wrap(~stim)

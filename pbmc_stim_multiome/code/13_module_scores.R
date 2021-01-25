library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(viridis)
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

pbmc_lll_t <- subset(pbmc_lll, predicted.celltype.l1 %in% c("other T", "CD4 T", "CD8 T"))
dim(pbmc_lll_t)

write.table(pbmc_lll@meta.data[,c("stim","protein_MS1", "rna_MS1","peaks.MS", "predicted.celltype.l1")],
            file = "../output/LLL_module_scores.csv", row.names = TRUE, col.names = TRUE, sep = ",", quote = FALSE)

FeaturePlot(pbmc_lll_t, features = c("peaks.MS","rna_MS1", "protein_MS1"),
                       min.cutoff = "q10", max.cutoff = "q90", split.by = "stim",
                       reduction =  'wnn.3.umap',  pt.size = 0.1,ncol = 3,
                       by.col = FALSE) &
  scale_color_viridis() &
  theme_void() & ggtitle("") &
  theme(legend.position = "none") &
  theme(plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"))

cowplot::ggsave2(p_MS_go, file = "../plots/umap_3wnn_tcellactivation_modulescores.png", width = 9, height = 6)

values_df <- data.frame(data.matrix(pbmc_lll_t@meta.data[pbmc_lll_t@meta.data$stim == "Stim",c("protein_MS1", "rna_MS1","peaks.MS")]))
library(gg3D)
values_df$peaksMS <- ifelse(values_df$peaks.MS > 20, 20,  ifelse(values_df$peaks.MS < -10, -10, values_df$peaks.MS))
p1 <- ggplot(values_df, aes(x=rna_MS1, y=protein_MS1, color=peaksMS)) + 
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "none") +
  labs(x = "RNA Activation Module", y = "Protein Activation Module") +
  geom_point(size = 0.2) +
  scale_color_viridis(limits = c(-10, 20))
cowplot::ggsave2(p1, file = "../plots/mod3_view_score_go.pdf", width = 1.6, height = 1.6)

ggplot(pbmc_lll_t@meta.data, aes(x = rna_MS1, y = protein_MS1, color = peaks.MS)) +
  geom_point() + facet_wrap(~stim)

library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BuenColors)
library(Signac)

pbmc_lll <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/output/pbmc_LLL_processed.rds")

# Process annotations
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
names(annotations) <- annotations$gene_name
pbmc_lll@assays[["peaks"]] <- RegionStats(pbmc_lll@assays[["peaks"]], genome = BSgenome.Hsapiens.UCSC.hg38)

gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annotations)

map <- fread("../../pbmc_stimulation_asapseq_citeseq/data/marker_gene_mapping.tsv")
map <- map[(map$Marker_name %in% rownames(pbmc_lll@assays$ADT@data)) & (!is.na(map[[2]])),]
data <- pbmc_lll@assays$ADT@data[map[[1]] ,]
rownames(data) <- map[[2]]
pbmc_lll@assays[["ADTforP2G"]] <- CreateAssayObject(data = data)

rna_peak_gene <- LinkPeaks(
  pbmc_lll, "peaks", "RNA",
  expression.slot = "data",
  distance = 5e+05,
  min.cells = 10,
  method = "pearson",
  genes.use = map[[2]],
  n_sample = 200,
  pvalue_cutoff = 1,
  gene.coords = gene.coords,
  score_cutoff = 0,
  verbose = TRUE
)@assays$peaks@links %>% data.frame()

adt_peak_gene <- LinkPeaks(
  pbmc_lll, "peaks", "ADTforP2G",
  expression.slot = "data",
  distance = 5e+05,
  min.cells = 10,
  method = "pearson",
  genes.use = map[[2]],
  n_sample = 200,
  pvalue_cutoff = 1,
  gene.coords = gene.coords,
  score_cutoff = 0,
  verbose = TRUE
)@assays$peaks@links %>% data.frame()

mdf <- merge(rna_peak_gene, adt_peak_gene, by = c("peak", "gene"))

p1 <- ggplot(mdf, aes(x = zscore.x, y = zscore.y, label = gene)) +
  geom_point(size = 0.2, alpha = 0.5)+ pretty_plot(fontsize = 7) + L_border() +
  labs(x = "Peak to RNA z-score ", y= "Zscore Peak to Protein")
cowplot::ggsave2(p1, file = "../plots/peak_to_gene_zscores.pdf", width = 1.8, height = 1.8)

library(Seurat)
library(data.table)
library(Signac)
library(dplyr)
library(Matrix)
library(harmony)
library(viridis)
library(dplyr)
library(BuenColors)
library(EnsDb.Hsapiens.v86)

# Import protein
import_kite_counts <- function(lib_id, barcode_idx = ""){
  mtx <- fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.barcodes.txt"), header = FALSE)[[1]], "-", barcode_idx)
  colnames(matx) <- paste0(fread(paste0("../data/protein_data/",lib_id,"_kite/featurecounts/featurecounts.genes.txt"), header = FALSE)[[1]])
  return(t(matx))
}

all_protein <- cbind(import_kite_counts("LLL_ctrl", "1"),import_kite_counts("LLL_stim", "2"), import_kite_counts("DIG_ctrl", "3"), import_kite_counts("DIG_stim", "4"))

# Import RNAseq
# Import scRNA-seq
process_scRNA  <- function(raw, new_idx){
  
  # Remove crazy high and low expressors
  n_feature_rna <- colSums(raw > 0)
  n_total_rna <- colSums(raw)
  pct_mito <- colSums(raw[grepl("^MT", rownames(raw)), ])/n_total_rna * 100
  qc_cells <- colnames(raw)[pct_mito < 30 & n_total_rna > 1000 & n_feature_rna > 500]
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), qc_cells]
  colnames(raw) <- gsub("-1", paste0("-", new_idx), colnames(raw))
  raw
}

# Import the RNA
lll_ctrl_rna <- process_scRNA(Read10X_h5("../data/LLL_CTRL_feature_bc_matrix.h5")$`Gene Expression`, "1")
lll_stim_rna <- process_scRNA(Read10X_h5("../data/LLL_STIM_feature_bc_matrix.h5")$`Gene Expression`, "2")
dig_ctrl_rna <- process_scRNA(Read10X_h5("../data/DIG_CTRL_feature_bc_matrix.h5")$`Gene Expression`, "3")
dig_stim_rna <- process_scRNA(Read10X_h5("../data/DIG_STIM_feature_bc_matrix.h5")$`Gene Expression`, "4")

# Import the meta data and create seurat object
import_qc <- function(file, new_idx){
  df <- fread(file) %>% dplyr::filter(is_cell == 1) %>% data.frame()
  df$barcode <- gsub("-1", paste0("-", new_idx), df$barcode)
  rownames(df) <-df$barcode
  df
}
tenx_qc <- rbind(
  import_qc("../data/metrics/LLL_CTRL_per_barcode_metrics.csv.gz","1"),
  import_qc("../data/metrics/LLL_STIM_per_barcode_metrics.csv.gz","2"),
  import_qc("../data/metrics/DIG_CTRL_per_barcode_metrics.csv.gz","3"),
  import_qc("../data/metrics/DIG_STIM_per_barcode_metrics.csv.gz","4")
)

full_meta <- merge(tenx_qc, data.frame(barcode = colnames(all_protein), 
                                       CD4adt = all_protein["CD4-1",],
                                       CD8adt = all_protein["CD8",],
                                       totalADT = colSums(all_protein), 
                                       totalCTRLadt = colSums(all_protein[grepl("Ctrl", rownames(all_protein)), ])), by = "barcode")
rownames(full_meta) <- full_meta$barcode
dim(full_meta)

rownames(full_meta) <- full_meta$barcode
dim(full_meta)
full_meta$channel <- substr(rownames(full_meta), 18, 18)
table(full_meta$channel[full_meta$totalCTRLadt < 10] )
table(full_meta$channel[full_meta$totalCTRLadt < 10] )
full_meta %>% dplyr::filter(channel == "3") %>% mutate(rat=atac_peak_region_fragments/atac_fragments*100) %>% pull(rat )%>% length()

# Do some visualization to determine thresholds
ggplot(full_meta, aes(x = totalCTRLadt, y = totalADT)) +
  geom_point() + scale_y_log10() + scale_x_log10()  +
  facet_wrap(~channel)
sum(full_meta$totalADT > 1e4 | full_meta$totalCTRLadt > 9)
full_meta$pct_in_peaks <- full_meta$atac_peak_region_fragments / full_meta$atac_fragments *100

pbmc_all <- CreateSeuratObject(counts = cbind(lll_ctrl_rna,lll_stim_rna, dig_ctrl_rna, dig_stim_rna), 
                               meta.data = full_meta)
qplot(pbmc_all$pct_in_peaks)

# Prospectively subset on attributes
pbmc_all <- subset(pbmc_all, subset = pct_in_peaks > 50 & totalCTRLadt < 10 & totalADT > 100 & nCount_RNA < (10^4.5) & !(CD8adt > 30 & full_meta$CD4adt > 100))
pbmc_all@meta.data$stim <- ifelse(substr(colnames(pbmc_all), 18, 18) %in% c(1,3), "Control", "Stim")
pbmc_all@meta.data$assay <- ifelse(substr(colnames(pbmc_all), 18, 18) %in% c(1,2), "LLL", "DIG")
pbmc_all@meta.data$stim_assay <- paste0(pbmc_all@meta.data$stim, "_", pbmc_all@meta.data$assay)
table(pbmc_all@meta.data$stim_assay)

# Import ATAC-seq
peaks <- Read10X_h5("../../../asap_large_data_files/multiome_pbmc_stim/input/filtered_peak_bc_matrix.h5")
lll_peaks <- peaks[,colnames(pbmc_all)]

chrom_assay <- CreateChromatinAssay(
  counts = lll_peaks,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '../../../asap_large_data_files/multiome_pbmc_stim/input/fragments.tsv.gz',
  min.cells = 0,
  min.features = 0
)

pbmc_all[["peaks"]] <- chrom_assay
pbmc_all[["ADT"]] <- CreateAssayObject(all_protein[,colnames(pbmc_all)])

# RNA Seurat stuff
DefaultAssay(pbmc_all) <- "RNA"
pbmc_all <- pbmc_all  %>% 
  NormalizeData() %>% ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA(verbose = FALSE, assay = "RNA", reduction.name = "pca") %>%
  RunHarmony( group.by.vars = c('stim',"assay"), reduction = 'pca', assay.use = 'RNA',project.dim = FALSE,  reduction.save = "harmony_RNA")

# LSI dim reduction
DefaultAssay(pbmc_all) <- "peaks"
pbmc_all <- RunTFIDF(pbmc_all) %>% 
  FindTopFeatures( min.cutoff = 'q25') %>%
  RunSVD() %>%
  RunHarmony( group.by.vars = c('stim',"assay"), reduction = 'lsi', assay.use = 'peaks',project.dim = FALSE,  reduction.save = "harmony_Peaks")

# Do it for ADT
DefaultAssay(pbmc_all) <- "ADT"
pbmc_all <- pbmc_all  %>% 
  NormalizeData(assay = "ADT", normalization.method = "CLR", margin = 2) %>%
  ScaleData(assay = "ADT", do.scale = FALSE) %>%
  FindVariableFeatures(assay = "ADT") %>% 
  RunPCA(verbose = FALSE, assay = "ADT", reduction.name = 'apca') %>%
  RunHarmony( group.by.vars = c('stim',"assay"), reduction = 'apca', assay.use = 'ADT',project.dim = FALSE,  reduction.save = "harmony_ADT")

# Now run multimodal neighbors and embedding
pbmc_all <- FindMultiModalNeighbors(object = pbmc_all,
                                    reduction.list = list("harmony_RNA", "harmony_Peaks", "harmony_ADT"),
                                    dims.list = list(1:30, 2:30, 1:30))
pbmc_all <- RunUMAP(pbmc_all, nn.name = "weighted.nn", reduction.name = "wnn.3.umap", reduction.key = "Uw3_" )
pbmc_all <- FindClusters(pbmc_all, graph.name = "wsnn", algorithm = 3, resolution = 0.2, verbose = FALSE) # 
pbmc_all <- FindClusters(pbmc_all, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

# Visualize
p1 <- DimPlot(pbmc_all, reduction = 'wnn.3.umap', label = FALSE, pt.size = 0.01, split.by = "stim_assay" , group.by = "assay") & 
  scale_color_manual(values = c("dodgerblue3", "firebrick")) & theme_void() &
  theme(legend.position = "none")
cowplot::ggsave2(p1, file = "../plots/dig_lll_out_all4.png", width = 14, height = 3.6, dpi = 500)

# Visualize some stuff
DefaultAssay(pbmc_all) <- "ADT"
p_feature <- FeaturePlot(pbmc_all, features = c( "CD138-1(Syndecan-1)"), split.by = "stim",
                         reduction =  'wnn.3.umap',  min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.01) &
  theme_void() & scale_color_gradientn(colors = jdb_palette("solar_extra")) &
  theme(legend.position = "none")
cowplot::ggsave2(p_feature, file = "../plots/cd138_expression_LLLdig.pdf", 
                 width = 3.5, height = 1.6)

p_embedding <- DimPlot(pbmc_all, reduction = 'wnn.3.umap', label = TRUE, repel = TRUE, label.size = 2, 
                       pt.size = 0.01,
                       split.by = "stim")  + theme_void() + theme(legend.position = "none")  +
  scale_color_manual(values = jdb_palette("corona"))

cowplot::ggsave2(p_embedding, file = "../plots/stim_split_embedding.pdf", 
                 width = 4.5, height = 1.8)


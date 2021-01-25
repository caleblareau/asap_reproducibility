library(Seurat)
library(BuenColors)
library(data.table)
library(dplyr)
# the 10x hdf5 file contains both data types. 
import_pctMT_exon_RNAseq <- function(file,qc_file,  condition, what){
  inputdata.10x <- Read10X_h5(file)
  pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`)
  initial_df <- merge(data.frame(
    barcode = colnames(pbmc),
    what, 
    condition, 
    count_mt = PercentageFeatureSet(pbmc, pattern = "^MT-")[[1]] * pbmc@meta.data$nCount_RNA /100),
    fread(qc_file) %>% filter(is_cell == 1), by = "barcode")
  
  initial_df[,c("count_mt", "gex_exonic_umis", "gex_intronic_umis", "what", "condition")] %>%
    mutate(exon_no_mt = gex_exonic_umis-count_mt) %>% 
    mutate(pct_exon = gex_exonic_umis/(gex_exonic_umis + gex_intronic_umis)*100) %>%
    mutate(pct_exon_no_mt = exon_no_mt/(exon_no_mt + gex_intronic_umis)*100) 
  
}

all_mtRNA_df <- rbind(
  import_pctMT_exon_RNAseq("../../../asap_large_data_files/multiome_pbmc_stim/input/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5",
                      '../data/metrics/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv.gz',"Ctrl", "Multi"),
  import_pctMT_exon_RNAseq("../data/DIG_CTRL_feature_bc_matrix.h5", 
                      "../data/metrics/DIG_CTRL_per_barcode_metrics.csv.gz", "Ctrl", "DIG"),
  import_pctMT_exon_RNAseq("../data/DIG_STIM_feature_bc_matrix.h5", 
                      "../data/metrics/DIG_STIM_per_barcode_metrics.csv.gz", "Stim", "DIG"),
  import_pctMT_exon_RNAseq("../data/LLL_CTRL_feature_bc_matrix.h5", 
                      "../data/metrics/LLL_CTRL_per_barcode_metrics.csv.gz", "Ctrl", "LLL"),
  import_pctMT_exon_RNAseq("../data/LLL_STIM_feature_bc_matrix.h5", 
                      "../data/metrics/LLL_STIM_per_barcode_metrics.csv.gz", "Stim", "LLL")
)

propRNA1 <- ggplot(all_mtRNA_df, aes(x = condition, y = (pct_exon), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mapped to exons", color = "") +
  theme(legend.position = "none")

propRNA2 <- ggplot(all_mtRNA_df, aes(x = condition, y = (pct_exon_no_mt), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs mappted exons - no mito", color = "") +
  theme(legend.position = "none")


cowplot::ggsave2(cowplot::plot_grid(
  propRNA1, propRNA2, nrow =1
), file = "../plots/RNA_features_pct_QC.pdf", width = 3.3, height = 1.5)

all_mtRNA_df %>% group_by(condition, what) %>% summarize(median(pct_exon), median(pct_exon_no_mt))

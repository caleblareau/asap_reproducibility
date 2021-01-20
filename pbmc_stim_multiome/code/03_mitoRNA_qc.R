library(Seurat)
library(BuenColors)

# the 10x hdf5 file contains both data types. 
import_pctMT_RNAseq <- function(file, condition, what){
  inputdata.10x <- Read10X_h5(file)
  pbmc <- CreateSeuratObject(counts = inputdata.10x$`Gene Expression`)
  data.frame(
    what, 
    condition, 
    pct_mt = PercentageFeatureSet(pbmc, pattern = "^MT-")[[1]])
}

all_mtRNA_df <- rbind(
  import_pctMT_RNAseq("../../../asap_large_data_files/multiome_pbmc_stim/input/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5", "Ctrl", "Multi"),
  import_pctMT_RNAseq("../data/DIG_CTRL_feature_bc_matrix.h5", "Ctrl", "DIG"),
  import_pctMT_RNAseq("../data/DIG_STIM_feature_bc_matrix.h5", "Stim", "DIG"),
  import_pctMT_RNAseq("../data/LLL_CTRL_feature_bc_matrix.h5", "Ctrl", "LLL"),
  import_pctMT_RNAseq("../data/LLL_STIM_feature_bc_matrix.h5", "Stim", "LLL")
)

pMitoRNA <- ggplot(all_mtRNA_df, aes(x = condition, y = (pct_mt), color = what)) +
  geom_boxplot(outlier.shape = NA,  position = position_dodge(preserve = "single")) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,50)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%UMIs from mtRNA", color = "") +
  theme(legend.position = "none")
  
pMitoRNA


cowplot::ggsave2(cowplot::plot_grid(
  pMitoRNA
), file = "../plots/mitoRNApct_QC.pdf", width = 1.65, height = 1.5)

all_mtRNA_df %>% group_by(condition, what) %>% summarize(median(pct_mt))

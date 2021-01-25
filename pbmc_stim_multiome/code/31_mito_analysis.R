library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

source("../../global_functions/variant_calling.R")

# Expensive computation to ID the variants
if(FALSE){
  
  # Pull out the cells
  ctrl <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/input/LLL_CTRL_GRCh38-mtMask_atac_mgatk.rds") 
  stim <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/input/LLL_STIM_GRCh38-mtMask_atac_mgatk.rds")
  table(stim$depth > 20)
  table(ctrl$depth > 20)
  
  colnames(stim) <- gsub("-1", "-2", colnames(stim))
  both_mgatk <- call_mutations_mgatk(cbind(ctrl[,ctrl$depth > 20], stim[,stim$depth > 20]))
  
  # Make the standard mtscatac variant calling plot
  misc_df <- data.frame(rowData(both_mgatk))
  p1 <- ggplot(misc_df %>%  filter(n_cells_conf_detected >= 5 ), aes(x = strand_correlation, y = log10(vmr), color = log10(vmr) > -2 & strand_correlation > 0.65)) +
    geom_point(size = 0.4) + scale_color_manual(values = c("black", "firebrick")) +
    labs(color = "HQ", x = "Strand concordance", y = "log VMR") +
    pretty_plot(fontsize = 8) + L_border() +
    theme(legend.position = "bottom") + 
    geom_vline(xintercept = 0.65, linetype = 2) +
    geom_hline(yintercept = -2, linetype = 2) + theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "LLLboth_var_call_mtscatac.pdf", width = 1.55, height = 1.55)
  saveRDS(assays(both_mgatk)[["allele_frequency"]][variants,], file = "both_LLL_allele_freqs.rds")
}

af <- readRDS("../output/both_LLL_allele_freqs.rds")
set.seed(1)
stimvec <- substr(colnames(af), 18, 18) == 2
stimvec_permuted <- sample(stimvec)
af_change_df <- data.frame(
  variant = rownames(af),
  observed = log2(rowMeans(af[,stimvec])/rowMeans(af[,!stimvec])),
  permuted = log2(rowMeans(af[,stimvec_permuted])/rowMeans(af[,!stimvec_permuted]))
)

ks.test(af_change_df$observed, af_change_df$permuted)
density_df <- af_change_df %>% reshape2::melt()
p2 <- ggplot(density_df, aes(x = value, fill = variable)) + 
  geom_density(alpha = 0.5) +
  scale_fill_manual(values =c ("firebrick", "grey")) +
  scale_y_continuous(expand = c(0,0)) +
  pretty_plot(fontsize = 8) + L_border() +
  labs(x = "log2 FC stim / control", y = "empirical density") +
  theme(legend.position = "none")

cowplot::ggsave2(p2, file = "../plots/allele_frequency_shift_stim.pdf", width = 1.9, height = 1.6)

# Look at the global gain of mtDNA
ctrl <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/input/LLL_CTRL_GRCh38-mtMask_atac_mgatk.rds") 
stim <- readRDS("../../../asap_large_data_files/multiome_pbmc_stim/input/LLL_STIM_GRCh38-mtMask_atac_mgatk.rds")


posdf <- data.frame(
  pos = 1:16569,
  ctrl = rowMeans(assays(ctrl)[["coverage"]]),
  stim = rowMeans(assays(stim)[["coverage"]])
)
pcoverage  <- reshape2::melt(posdf , id.vars = "pos") %>%
  dplyr::filter(value > 25)%>%  # one outlier point
  ggplot(aes(x = pos, y = value, color = variable)) +
  geom_line() +
  pretty_plot(fontsize = 7) + L_border() +
  scale_color_manual(values = c("dodgerblue3", "firebrick"))  +
  theme(legend.position = "none") + labs(x = "Position along mtDNA genome", y = "Mean coverage of mtDNA per cell")
cowplot::ggsave2(pcoverage, file = "../plots/Stim_mtDNAcoverage.pdf", width = 2.5, height = 1.5)

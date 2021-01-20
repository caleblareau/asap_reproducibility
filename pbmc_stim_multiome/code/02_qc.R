library(dplyr)
library(data.table)
library(stringr)
library(BuenColors)

source('../../global_functions/estimateLibraryComplexity.R')

process_qc <- function(dt, what ){
  dt %>% filter(is_cell == 1) %>% arrange(desc(atac_fragments)) %>% 
    mutate(pct_mito = atac_mitochondrial_reads/ atac_raw_reads) %>%
    mutate(assay = what)  %>%
    mutate(atac_total = atac_dup_reads + atac_fragments) -> odf
  odf$atac_complexity <- sapply(1:dim(odf)[1], function(i){
   estimateLibrarySize(odf$atac_total[i] ,odf$atac_fragments[i])
  })
  odf[,c("pct_mito", "atac_complexity", "gex_genes_count","assay")]
}

dig <- process_qc(fread('../data/metrics/DIG_CTRL_per_barcode_metrics.csv.gz'), "DIG")
multi <- process_qc(fread('../data/metrics/pbmc_granulocyte_sorted_10k_per_barcode_metrics.csv.gz'), "Multi")
LLL <-  process_qc(fread('../data/metrics/LLL_CTRL_per_barcode_metrics.csv.gz'), "LLL")

all_df_qc <- rbind(dig, multi, LLL)

pMito <- ggplot(all_df_qc, aes(x = assay, y =pct_mito*100, color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%mtDNA (ATAC)", color = "") +
  theme(legend.position = "none")

pATAC <- ggplot(all_df_qc, aes(x = assay, y = log10(atac_complexity), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 ATAC Complexity", color = "") +
  theme(legend.position = "none")

pRNA <- ggplot(all_df_qc, aes(x = assay, y = (gex_genes_count), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 genes detected (RNA)", color = "") +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(
  pATAC, pMito, pRNA, nrow = 1
), file = "../plots/all_ctrl_QC.pdf", width = 4, height = 1.5)

all_df_qc %>% group_by(assay) %>%
  summarize(mito = median(pct_mito*100), median(atac_complexity), 
            median(gex_genes_count), count = n())

digS <- process_qc(fread('../data/metrics/DIG_STIM_per_barcode_metrics.csv'), "DIG")
LLLS <-  process_qc(fread('../data/metrics/LLL_STIM_per_barcode_metrics.csv'), "LLL")

all_df_qc_stim <- rbind(digS, LLLS)

pMitoS <- ggplot(all_df_qc_stim, aes(x = assay, y =pct_mito*100, color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(0,60)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "%mtDNA (ATAC)", color = "") +
  theme(legend.position = "none")

pATACS <- ggplot(all_df_qc_stim, aes(x = assay, y = log10(atac_complexity), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5.5)) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 ATAC Complexity", color = "") +
  theme(legend.position = "none")

pRNAS <- ggplot(all_df_qc_stim, aes(x = assay, y = (gex_genes_count), color = assay)) +
  geom_boxplot(outlier.shape = NA) + 
  pretty_plot(fontsize = 7) + L_border() + scale_y_log10(limits = c(10^2, 10^4), breaks = c(100, 500, 1000, 5000, 10000) ) +
  scale_color_manual(values = c("dodgerblue3", "firebrick", "darkgrey")) +
  labs(x = "Condition", y = "log10 genes detected (RNA)", color = "") +
  theme(legend.position = "none")

cowplot::ggsave2(cowplot::plot_grid(
  pATACS, pMitoS, pRNAS, nrow = 1
), file = "../plots/all_stim_QC.pdf", width = 3.5, height = 1.5)

all_df_qc_stim %>% group_by(assay) %>%
  summarize(mito = median(pct_mito*100), median(atac_complexity), 
            median(gex_genes_count), count = n())


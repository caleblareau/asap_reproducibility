library(BuenColors)
library(dplyr)
source("../../global_functions/variant_enrichment.R")

af <- readRDS("../output/both_LLL_allele_freqs.rds")
ref_all <- fread("../../global_functions/references/chrM_refAllele.txt")

prop_df <- get_enrich_mutation_df(rownames(af), ref_all)

# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide (n = 106 mutations)", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "../plots/multiome_n106_mito_signature.pdf", width = 4, height = 2.4)

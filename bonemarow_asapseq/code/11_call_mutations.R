library(SummarizedExperiment)
library(dplyr)

source("../../global_functions/variant_calling.R")
source("../../global_functions/variant_enrichment.R")

# Import and filter cells
SE <- readRDS("../../../asap_large_data_files/bonemarrow_data/input/BM_asap_v12_hg38-mtMask_mgatk.rds")
barcodes <- fread("../data/barcodes/step3_ADThq.tsv", header = FALSE)[[1]]
SE_filt <- SE[,colnames(SE) %in% barcodes & colData(SE)$depth > 10]

# Export a fasta of the most likely genotype / position
if(FALSE){
  library(seqinr)
  get_letter_total <- function(letter) rowSums(assays(SE)[[paste0(letter, "_counts_fw")]]) + rowSums(assays(SE)[[paste0(letter, "_counts_rev")]])
  mat <- data.matrix(data.frame(
    xA = get_letter_total("A"),
    xC = get_letter_total("C"),
    xG = get_letter_total("G"),
    xT = get_letter_total("T")
  ))
  chars <- c("A", "C", "G", "T")[(max.col(mat))]
  haplotype <- paste(chars, collapse = "")
  write.fasta(sequences = haplotype, names = "chrM", file.out = "../output/chrM_haplotype_donor_16260C.fasta")
  
  # Look at the T
  chars_T <- chars
  chars_T[16260] <- "T"
  write.fasta(sequences =  paste(chars_T, collapse = ""), names = "chrM", file.out = "../output/chrM_haplotype_donor_16260T.fasta")
  
}

# Now call mutations
mSE <- call_mutations_mgatk(SE_filt)
muts_df <- rowData(mSE)
boo <- muts_df$n_cells_conf_detected >= 3 & muts_df$strand_correlation > 0.65 & log10(muts_df$vmr) > -2 & muts_df$mean_coverage >= 10
vars_clone <- data.frame(muts_df) %>%  filter(boo) 
vars_clone %>% arrange(desc(mean))

dim(vars_clone)
saveRDS(mSE[boo,], file = "../output/mitoMutations_cov10_marrow.rds")


ref_all <- fread("../../global_functions/references/chrM_refAllele.txt")

prop_df <- get_enrich_mutation_df(vars_clone$variant, ref_all)

# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide (n = 99 mutations)", y = "Substitution Rate (Expected / Observed)")
cowplot::ggsave2(p1, file = "../plots/bonemarrow_n99_mito_signature.pdf", width = 4, height = 2.4)

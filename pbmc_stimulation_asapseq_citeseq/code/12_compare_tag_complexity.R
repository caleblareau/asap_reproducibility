library(dplyr)
library(data.table)
library(stringr)

pass_bc_df <- data.frame(fread("../output/HQ_barcodes_4exps.tsv"))

# Assign the bus set to the corresponding 
process_tag_reference <- function(ref_file){
  
  ref_df <- fread(ref_file, header = FALSE)
  ref_df2 <- data.frame(
    seq = ref_df[[1]][c(FALSE,TRUE)],
    name = gsub(">", "", ref_df[[1]][c(TRUE,FALSE)])
  )
  ref_df2$idx <- as.character(1:dim(ref_df2)[1])
  ref_df2$name_short <- str_split_fixed(ref_df2$name, "-", 2)[,1]
  
  # Some of these Abs have an extra -# that is necessary
  doubles <- data.frame(table(ref_df2$name_short)) %>% filter(Freq > 50) %>% pull(Var1) %>% as.character()
  ref_df2$extra_idx <- str_split_fixed(ref_df2$name, "-", 3)[,2]
  ref_df2$use <- ifelse(ref_df2$name_short %in% doubles, paste0(ref_df2$name_short,"_", ref_df2$extra_idx),ref_df2$name_short)
  
  vec <- as.character(ref_df2$use); names(vec) <- as.character(ref_df2$idx)
  vec
}

# Import PBMC references
pbmc_vec <- process_tag_reference("../data/references/PBMC_FeaturesMismatch.fa")

source('../../global_functions/estimateLibraryComplexity.R')

# Function to process the per-barcode summary statistics
estimate_tag_rates <- function(bus_file, bcs, vec){
  bcs <- substr(bcs, 1, 16)
  dff <- fread(bus_file, col.names = c("barcode", "umi", "set", "idx")) %>%
    mutate(Ab = vec[as.character(set)], umi_ab = paste0(Ab,"_",umi)) %>%
    filter(barcode %in% bcs) %>% 
    group_by(barcode) %>%
    summarise(count = n_distinct(umi_ab), n_total = n()) %>%
    mutate(n_duplicated = n_total - count) 
  dff$complexity <- sapply(1:(dim(dff)[1]), function(i){
    estimateLibrarySize(dff$n_total[i],dff$count[i])
  })
  dff
}

# Pull HQ barcodes from embedding 
Ac_bcs <- pass_bc_df %>% filter(assay == "ATAC_control") %>% pull(barcode) %>% as.character()
As_bcs <- pass_bc_df %>% filter(assay == "ATAC_stim") %>% pull(barcode) %>% as.character()
Cc_bcs <- pass_bc_df %>% filter(assay == "RNA_control") %>% pull(barcode) %>% as.character()
Cs_bcs <- pass_bc_df %>% filter(assay == "RNA_stim") %>% pull(barcode) %>% as.character()

# Compute per-barcode statistics
ASAP_ctrl <- estimate_tag_rates("../../../asap_large_data_files/pbmc_stim_data/input/bus_files/ASAP_ctrl_featurecounts.tsv", Ac_bcs, pbmc_vec)
ASAP_stim <- estimate_tag_rates("../../../asap_large_data_files/pbmc_stim_data/input/bus_files/ASAP_stim_featurecounts.tsv", As_bcs, pbmc_vec)
CITE_ctrl <- estimate_tag_rates("../../../asap_large_data_files/pbmc_stim_data/input/bus_files/CITE_ctrl_featurecounts.tsv", Cc_bcs, pbmc_vec)
CITE_stim <- estimate_tag_rates("../../../asap_large_data_files/pbmc_stim_data/input/bus_files/CITE_stim_featurecounts.tsv", Cs_bcs, pbmc_vec)

save(ASAP_ctrl, ASAP_stim, CITE_ctrl, CITE_stim, file = "../output/complexity_matrices.rda")


if(FALSE){
  
  load("../output/complexity_matrices.rda")

  compute_duplicate_percentage <- function(df){
    df %>% filter(complexity < 100000) %>%
      mutate(percent_duplicates = n_duplicated/n_total*100) %>%
      pull(percent_duplicates) %>% summary
  }
  compute_duplicate_percentage(ASAP_ctrl)
  compute_duplicate_percentage(CITE_ctrl)
  compute_duplicate_percentage(ASAP_stim)
  compute_duplicate_percentage(CITE_stim)
  
  ASAP_ctrl$assay <- "ASAP";   ASAP_ctrl$condition <- "control" 
  CITE_ctrl$assay <- "CITE";   CITE_ctrl$condition <- "control" 
  ASAP_stim$assay <- "ASAP";   ASAP_stim$condition <- "stim" 
  CITE_stim$assay <- "CITE";   CITE_stim$condition <- "stim" 
  
  all_df <- rbind(ASAP_ctrl, CITE_ctrl, ASAP_stim, CITE_stim)
  p1 <- ggplot(all_df, aes(x = condition, y = log10(complexity), color = assay)) +
    geom_boxplot(outlier.shape = NA) + 
    pretty_plot(fontsize = 8) + L_border() + scale_y_continuous(limits = c(2.5, 4.5)) +
    scale_color_manual(values = c("firebrick", "dodgerblue3")) +
    labs(x = "Condition", y = "log10 Tag Library Complexity", color = "") +
    theme(legend.position = "bottom")
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay.pdf", width = 3, height = 2.4)
  
  all_df %>% group_by(assay, condition) %>% summarize(median(complexity))
}



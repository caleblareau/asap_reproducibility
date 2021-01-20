library(dplyr)
library(data.table)
library(stringr)


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
pbmc_vec <- process_tag_reference("../data/FeaturesMismatch.fa")

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
Dc <- fread("../data/metrics/DIG_CTRL_per_barcode_metrics.csv") %>% filter(is_cell == 1) %>% pull(barcode) %>% as.character()
Ds <- fread("../data/metrics/DIG_STIM_per_barcode_metrics.csv") %>% filter(is_cell == 1) %>% pull(barcode) %>% as.character()
Lc <- fread("../data/metrics/LLL_CTRL_per_barcode_metrics.csv") %>% filter(is_cell == 1) %>% pull(barcode) %>% as.character()
Ls <- fread("../data/metrics/LLL_STIM_per_barcode_metrics.csv") %>% filter(is_cell == 1) %>% pull(barcode) %>% as.character()

# Compute per-barcode statistics
DIG_ctrl <- estimate_tag_rates("../DIG_ctrl_mo.bustools.tsv", Dc, pbmc_vec)
DIG_stim <- estimate_tag_rates("../DIG_stim_mo.bustools.tsv", Ds, pbmc_vec)
LLL_ctrl <- estimate_tag_rates("../LLL_ctrl_mo.bustools.tsv", Lc, pbmc_vec)
LLL_stim <- estimate_tag_rates("../LLL_stim_mo.bustools.tsv", Ls, pbmc_vec)

save(DIG_ctrl, DIG_stim, LLL_ctrl, LLL_stim, file = "../output/complexity_matrices.rda")


if(FALSE){
  
  load("../output/complexity_matrices.rda")

  compute_duplicate_percentage <- function(df){
    df %>% filter(complexity < 100000) %>%
      mutate(percent_duplicates = n_duplicated/n_total*100) %>%
      pull(percent_duplicates) %>% summary
  }
  compute_duplicate_percentage(DIG_ctrl)
  compute_duplicate_percentage(DIG_stim)
  compute_duplicate_percentage(LLL_ctrl)
  compute_duplicate_percentage(LLL_stim)
  
  DIG_ctrl$assay <- "DIG";   DIG_ctrl$condition <- "control" 
  DIG_stim$assay <- "DIG";   DIG_stim$condition <- "stim" 
  LLL_ctrl$assay <- "LLL";   LLL_ctrl$condition <- "control" 
  LLL_stim$assay <- "LLL";   LLL_stim$condition <- "stim" 
  
  all_df <- rbind(DIG_ctrl, DIG_stim, LLL_ctrl, LLL_stim)
  p1 <- ggplot(rbind(DIG_ctrl, LLL_ctrl), aes(x = assay, y = log10(complexity), color = assay)) +
    geom_boxplot(outlier.shape = NA) + 
    pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5)) +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "Condition", y = "log10 Tag Complexity", color = "") +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay_MO_ctrl.pdf", width = 1, height = 1.5)
  
  p1 <- ggplot(rbind(DIG_stim, LLL_stim), aes(x = assay, y = log10(complexity), color = assay)) +
    geom_boxplot(outlier.shape = NA) + 
    pretty_plot(fontsize = 7) + L_border() + scale_y_continuous(limits = c(2, 5)) +
    scale_color_manual(values = c("dodgerblue3", "firebrick")) +
    labs(x = "Condition", y = "log10 Tag Complexity", color = "") +
    theme(legend.position = "none")
  cowplot::ggsave2(p1, file = "../plots/tag_complexity_assay_MO_stim.pdf", width = 1, height = 1.5)
  
  
  all_df %>% group_by(assay, condition) %>% summarize(median(complexity), count = n())
}



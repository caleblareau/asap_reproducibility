library(Seurat)
library(viridis)
library(scales)
library(Matrix)
library(ComplexHeatmap)
library(BuenColors)
library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

# Import ATAC processed data
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Look at the ADT data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adt_simple_norm <- data.matrix(t(t(adt_ss)/colSums(adt_ss))*1000)
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")

mat <- data.matrix(adtbm@assays$ADT@scale.data)

# Compute mean to filter
m_ps <- mdf %>% dplyr::filter(!is.na(monocyte_PS)) %>% pull(barcode)
ps_mean_adt_count <- rowMeans(data.matrix(adt_simple_norm)[,m_ps])
names(ps_mean_adt_count) <- rownames(adt_simple_norm)
other_mean_adt_count <- rowMeans(data.matrix(adt_simple_norm)[,colnames(adt_simple_norm) %ni% m_ps])
names(other_mean_adt_count) <- rownames(adt_simple_norm)

# Process the gene scores and pseudotime axis
dat <- fread("../../multiple_datasets/data/marker_gene_mapping.tsv")
gs_mat <- readRDS("../../../asap_large_data_files/bonemarrow_data/output/signac_marrow_gene_scores.rds")

# Make a master data frame
master_df <- data.frame(
  dat, 
  adt.idx = match(dat$Marker_name,rownames(mat)),
  gs.idx = match(dat$Gene_symbol, rownames(gs_mat))
)
master_df$mean_adt <- ps_mean_adt_count[master_df$Marker_name]
master_df$mean_adt_other_cells <- other_mean_adt_count[master_df$Marker_name]

master_df <- master_df[complete.cases(master_df),] %>% mutate(ratio = mean_adt/mean_adt_other_cells)

# Order the cells and the gene / tags
ordered_cells <- mdf[!is.na(mdf$monocyte_PS),] %>% arrange((monocyte_PS))
gs_mat_ordered <- gs_mat[,ordered_cells$barcode]
adt_mat_ordered <- mat[,ordered_cells$barcode]

# Now group/order/smooth the cell states and make ArchR's nice heatmaps
# see: https://github.com/GreenleafLab/ArchR/blob/5f855d7b7ff3f57eb0d28f312c5ea5d373d27ebd/R/Trajectory.R
n_bins <- 100
groups_gs <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(monocyte_PS >= (idx-1)*multiplier & monocyte_PS < (idx*multiplier)) %>% pull(barcode)
  print(length(bcs))
  rs <- rowSums(gs_mat_ordered[,bcs, drop = FALSE], na.rm = TRUE)
  log2(rs/sum(rs) *10000 + 1) # basic normalization
})

groups_adt <- sapply(1:n_bins, function(idx){
  multiplier <- 100/n_bins
  bcs <- ordered_cells %>% dplyr::filter(monocyte_PS >= (idx-1)*multiplier & monocyte_PS < idx*multiplier) %>% pull(barcode)
  rowMeans(adt_mat_ordered[,bcs], na.rm = TRUE)
})

# Filter based on indices
groups_gs <- groups_gs[master_df$gs.idx,]
groups_adt <- groups_adt[master_df$adt.idx,]

smoothWindow = 11
smooth_groups_adt <- data.matrix((apply((groups_adt), 1, function(x) ArchR:::.centerRollMean((x), k = smoothWindow))))
smooth_groups_gs <- data.matrix((apply((groups_gs), 1, function(x) ArchR:::.centerRollMean((x), k = smoothWindow))))

smooth_groups_minmax_adt <- t(apply(smooth_groups_adt, 2, function(x)(x-min(x))/(max(x)-min(x))))
smooth_groups_minmax_gs <- t(apply(smooth_groups_gs, 2, function(x)(x-min(x))/(max(x)-min(x))))


time_df_50 <- data.frame(
  gs_50 = max.col(smooth_groups_minmax_gs >=0.5, ties.method = "first"),
  adt_50 = max.col(smooth_groups_minmax_adt >=0.5, ties.method = "first"),
  adt_max = max.col(smooth_groups_minmax_adt),
  gs_max = max.col(smooth_groups_minmax_gs),
  master_df
) %>% mutate(adt_after = adt_50 >= gs_50)

time_df_50 %>%
  filter(mean_adt > 1 & ratio > 1) %>%
  filter(adt_max > 25) %>%
  pull(adt_after) %>% table()

time_df_50 %>%
  filter(mean_adt > 1 & ratio > 1) %>%
  filter(adt_max > 25)  %>%

pS <- time_df_50 %>%
  filter(mean_adt > 1 & ratio > 1) %>%
  filter(adt_max > 25) %>% 
  ggplot(aes(x = gs_50, y = adt_50, color = gs_50 >= adt_50)) + geom_point(size = 0.5) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  pretty_plot(fontsize = 8) + L_border() +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme(legend.position = "none") +
  labs(x = "Pseudotime of 50% max gene score", y = "Pseudotime of 50% max protein expression") +
  scale_color_manual(values =c ("dodgerblue3", "firebrick"))
cowplot::ggsave2(pS, file = "../plots/adt_chromatin_timing_max50.pdf", width = 1.8, height = 1.8)  




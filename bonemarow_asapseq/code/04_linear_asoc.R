library(Seurat)
library(viridis)
library(scales)

# Import ATAC processed data
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Look at the ADT data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
mat <- data.matrix(adtbm@assays$ADT@scale.data)

corPS <- cor(mdf$myeloidPS, t(mat), use = "pairwise.complete", method = "spearman")
cor_pvalue <- sapply(1:dim(mat)[1], function(i){
  cor.test(mdf$myeloidPS, mat[i,], use = "pairwise.complete", method = "spearman")$p.value
})

data.frame(ADT = rownames(mat),
           corPS = corPS[1,],
           pvalue = cor_pvalue) %>% arrange(desc(corPS)) -> cor_df

write.table(cor_df, file = "../output/PS_cor_myeloid.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Make a lot of visuals
joint_df <- data.frame( mdf,t(mat)) 
joint_df <- joint_df[!is.na(joint_df$myeloidPS ), ]

lapply(20:261,function(i){
  ss_df <- joint_df[,c(17,18,16, i)]
  tag <- colnames(ss_df)[4]
  colnames(ss_df) <- c("UMAP1", "UMAP2", "PS", "tag")
  p1 <- ggplot(shuf(ss_df), aes(x = UMAP1, y = UMAP2, color = tag)) +
    geom_point() + scale_color_viridis(oob = squish, limits = c(0,3)) + 
    labs(color = tag) + pretty_plot() + L_border() + theme(legend.position = "bottom")
  
  p2 <- ggplot(shuf(ss_df), aes(x = PS, y = tag)) +
    geom_point(alpha = 1/5)  + geom_smooth() +
    labs(color = tag) + pretty_plot() + L_border() + theme(legend.position = "bottom") +
    labs(x = "Myeloid Pseudotime (Progenitor -> Monocyte)", y = paste0(tag, " CLR ADT expression"))
  
  ggsave(cowplot::plot_grid(p1, p2, nrow = 1), file = paste0("../plots/all_tags_MPS/18June_MyeloidPS_", tag, ".png"), 
         width = 8, height = 4)
  i
})


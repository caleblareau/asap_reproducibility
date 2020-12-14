library(ArchR)
library(BuenColors)
library(dplyr)
doubScores <- readRDS("../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38/ASAP_marrow_hg38/ASAP_marrow_hg38-Doublet-Summary.rds")
df <- data.frame(DS = doubScores$doubletResults$doubletScoreLSI,
                 DE = doubScores$doubletResults$doubletEnrichLSI, 
                 readRDS("../output/ArchR_main_metadata.rds"))

df$DE2 <- ifelse(df$DE > 3, 3, df$DE)
p1 <- ggplot(shuf(df) , aes(x = UMAP1, y = UMAP2, color = DE2)) +
  geom_point(size = 0.3) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  theme_void() + theme(legend.position = "none")
cowplot::ggsave2(p1, 
                 file = "../plots/umap_doublet.png", width = 1.8*4, height = 1.8*4, dpi = 300)

p2 <- qplot(df$DE, bins = 50) +
  labs(x = "Doublet Enrichment", y = "Count") +
  pretty_plot(fontsize = 8) + L_border()
cowplot::ggsave2(p2, file = "../plots/doublet_histogram.pdf", width = 1.8, height = 1.8)

# Look at associations between the tag and the doublet score

df$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(df))
counts <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,df$barcode]
mat <- log2(t(t(counts)/colSums(counts))*10000+1)


i <- 1
doub <- (df$DE2 >= 3)

lapply(1:dim(mat)[1], function(i){
  ddf <- data.frame(
    pvalue= wilcox.test(mat[i,doub],(mat[i,!doub]))$p.value,
    diff = log2(mean(mat[i,doub]) / mean((mat[i,!doub]))),
    factor = rownames(mat)[i], 
    what = "observed"
  )
 ddf
}) %>% rbindlist() %>% data.frame() %>% arrange(what) -> odf


lapply(1:dim(mat)[1], function(i){
  ct <- cor.test(df$DE, mat[i,])
  ddf <- data.frame(
    pvalue= ct$p.value,
    cor = ct$estimate,
    factor = rownames(mat)[i], 
    what = "observed"
  )
  set.seed(1)
  ct <- cor.test(sample(df$DE), mat[i,])
  permdf <- data.frame(
    pvalue= ct$p.value,
    cor = ct$estimate,
    factor = rownames(mat)[i], 
    what = "permuted"
  )
  
  rbind(permdf, ddf)
}) %>% rbindlist() %>% data.frame() %>% arrange(what, cor) -> odf


ggplot(odf, aes(y = -log10(pvalue), x = diff, color = what, label = factor)) +
  geom_text()

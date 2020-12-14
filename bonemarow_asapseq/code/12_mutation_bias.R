library(data.table)
library(dplyr)

archr_df <- readRDS(file = "../output/ArchR_main_metadata.rds")
archr_df$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
mSE <- readRDS("../output/mitoMutaitons_cov10_marrow.rds")
archr_cluster_vec <- archr_df$Clusters; names(archr_cluster_vec) <- as.character(archr_df$barcode) 

c_map<- c("Ery", "Ery", rep("Lym", 5), "Prog", # first 8
          rep("Lym", 5), rep("Myeloid", 8))
names(c_map) <- paste0("C", as.character(1:21))
clusters_ordered <- c_map[archr_cluster_vec[as.character(mSE@colData$sample)]]
var_meta_df <- data.frame(rowData(mSE))

# Make easy to interpet supplemental panels
data.frame(
  X13069GA = assays(mSE)[["allele_frequency"]]["13069G>A",],
  X13711GA = assays(mSE)[["allele_frequency"]]["13711G>A",],
  X16260CT = assays(mSE)[["allele_frequency"]]["16260C>T",],
  lineage = clusters_ordered
) %>%
  group_by(lineage) %>% 
  summarize(mX13069GA = mean(X13069GA > 0.05),
            mX13711GA = mean(X13711GA > 0.05),
            mX16260CT = mean(X16260CT > 0.5)) %>%
  reshape2::melt(id.vars = "lineage") -> for_hist_df

phist <- ggplot(for_hist_df, aes(x = lineage, y = value*100, fill = variable)) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  geom_bar(stat = "identity", color = "black") +
  pretty_plot(fontsize = 7) + labs(x = "Lineage", y = "% cells + for mutation") +
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("purple","green4", "blue"))
cowplot::ggsave2(phist, 
                 file = "../plots/bars_mtDNA.pdf", height = 1.8, width = 2.2)


# Do a per-variant anova
var_meta_df$kruskal_pvalue <- sapply(1:dim(var_meta_df)[1], function(i){
  var_df <- data.frame(
    cluster = clusters_ordered,
    variant = assays(mSE)[["allele_frequency"]][i,] 
  )
  res.kw <- kruskal.test(variant ~ cluster, data = var_df)
  (res.kw)[[3]]
})
var_meta_df$kruskal_pvalue_adj <- p.adjust(var_meta_df$kruskal_pvalue)
var_meta_df %>% arrange((kruskal_pvalue)) %>% head()

pMut_mean <- ggplot(var_meta_df, aes(x = -log10(kruskal_pvalue_adj), y = mean*100)) +
  geom_point(size = 0.5) + scale_y_log10(breaks = c(0.001, 0.01, 0.05, 0.1, 0.2, 0.5)*100) +
  labs(x = "-log10 p cluster bias", y = "Pseudobulk AF% (log10 scale)") +
  pretty_plot(fontsize = 8) + L_border()
#cowplot::ggsave2(pMut_mean, file = "../plots/mutation_bias_mean.pdf", width = 1.8, height= 1.8)


# Pull in UMAP and variants
var_df <- data.frame(barcode = colnames(mSE), clusters_ordered,
                     data.matrix(t(assays(mSE)[["allele_frequency"]])))
plot_df <- left_join(var_df, archr_df, by = "barcode")

pMain1 <- ggplot(shuf(plot_df), aes(x = UMAP1, y = UMAP2, color = X16260C.T > 0.5 )) + 
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "blue")) 

pMain2 <- ggplot(plot_df %>% arrange((X13069G.A)), aes(x = UMAP1, y = UMAP2, color = X13069G.A > 0.05 )) +
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "purple")) 

cowplot::ggsave2(cowplot::plot_grid(pMain1, pMain2, ncol =1), 
                 file = "../plots/two_mito_mutations_MF.png", height = 1.8*2, width = 0.9*2, dpi = 500)

pS1 <- ggplot(plot_df %>% arrange((X13711G.A)), aes(x = UMAP1, y = UMAP2, color = X13711G.A > 0.05 )) +
  geom_point(size = 0.2)  + theme_void() + ggtitle("") + theme(legend.position = "none") + scale_color_manual(values = c("lightgrey", "green4")) 

cowplot::ggsave2(pS1, 
                 file = "../plots/supplemental_mtDNA_plot.png", height = 1.8, width = 1.8, dpi = 500)

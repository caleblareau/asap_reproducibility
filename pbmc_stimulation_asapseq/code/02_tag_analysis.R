library(ArchR)
library(Seurat)

# Import ArchR Project file
path_arrow <- "../../../asap_large_data_files/pbmc_stim_data/output/archr_pbmc_stim/"
proj_file <- paste0(path_arrow, "ArchR_PBMCs_stim.rds")
proj <- readRDS(proj_file)

df <- data.frame(
  proj@embeddings$UMAP$df, 
  proj@cellColData$Clusters,
  proj@cellColData$Sample,
  proj@cellColData$DoubletEnrichment
)
colnames(df) <- c("UMAP1", "UMAP2", "cluster","sample", "DS")

set.seed(4)
umap_base <- ggplot(shuf(df), aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point(size = 0.1) + pretty_plot(fontsize = 7) + 
  facet_wrap(~sample) + scale_color_manual(values = sample(jdb_palette("corona")[1:15])) +
  theme_void() + theme(legend.position = "none") 

cowplot::ggsave2(umap_base, 
                 width = 6, height = 3, file = "../plots/umap_base_asap.png", dpi = 500)

# Import ADT counts
ctrl_adt <- readRDS("../output/adt_mat/ASAP_ctrl.rds")
stim_adt <- readRDS("../output/adt_mat/ASAP_stim.rds")
colnames(ctrl_adt) <- paste0("Control#", colnames(ctrl_adt), "-1")
colnames(stim_adt) <- paste0("Stim#", colnames(stim_adt), "-1")
adt <- cbind(ctrl_adt, stim_adt)
adt_ss <- adt[,rownames(df)]

# Process the ADT in Seurat
pbmc <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
pbmc <- NormalizeData(pbmc, assay = "ADT", normalization.method = "CLR")
pbmc <- ScaleData(pbmc, assay="ADT")
mat <- data.matrix(pbmc@assays$ADT@scale.data)

# Now impute the values
imputed_ADT_mat <- imputeMatrix(mat = mat, 
                                imputeWeights = getImputeWeights(proj),
                                logFile = "ArchRLogs/ArchR-addImputeWeights-RAND.log")

df_clr_scale <- data.frame(
  df, t(mat)
) 
df_archr_smooth <- data.frame(
  df, t(imputed_ADT_mat)
) 

saveRDS(df_clr_scale, file = "../output/adt_mat/ASAP_embedding_CLRadt.rds")
saveRDS(df_archr_smooth, file = "../output/adt_mat/ASAP_embedding_smoothadt.rds")


library(viridis)
library(scales)


ggplot(shuf(df_clr_scale), aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  geom_point() + scale_color_manual(values = sample(jdb_palette("corona")[1:15])) +
  facet_wrap(~sample)


(cor(df_clr_scale[,"CD158"],data.matrix(df_clr_scale[,c(-1,-2,-3,-4,-5)])))

ggplot(shuf(df_clr_scale), aes(x = UMAP1, y = UMAP2, color = CD158)) + 
  geom_point() + scale_color_viridis(oob = squish, limits = c(0,3)) + # 
  facet_wrap(~sample) + pretty_plot() 

shuf(df_clr_scale) %>%
  filter(cluster %in% c("C8", "C9")) %>%
  ggplot(aes(x = CD56.NCAM., y = CD2, color = sample)) +
  geom_point()

shuf(df_clr_scale) %>%
  filter(cluster %in% c("C8", "C9")) %>%
  ggplot(aes(x = CD56.NCAM., y = CD56.NCAM.Recombinant, color = sample)) +
  geom_point()


# Look at variability
var_df <- df_clr_scale[,c(-1,-2,-5)] %>%
  group_by(cluster, sample) %>%
  summarize_all(.funs = var) %>%
  reshape2::melt(id.vars = c('cluster', "sample")) %>%
  reshape2::dcast(cluster + variable ~ sample, value.var = "value")

ggplot(var_df, aes(x = Control, y = Stim, label = variable)) + 
  geom_text() + facet_wrap(~cluster, scales = "free")

nk_df <- data.frame(
  tag = rownames(mat),
  C8 = rowMeans(mat[,df$cluster == "C8"]),
  C9 = rowMeans(mat[,df$cluster == "C9"])
)
nk_df$diff <- nk_df$C8-nk_df$C9      
nk_df %>% arrange(diff) %>% tail()


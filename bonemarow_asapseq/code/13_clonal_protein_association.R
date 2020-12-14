library(data.table)
library(dplyr)
library(Seurat)
library(stringr)
library(SummarizedExperiment)
library(BuenColors)

"%ni%" <- Negate("%in%")
perc.rank <- function(x) trunc(rank(x))/length(x)

# Function to make cluster assignments from mtDNA heteroplasmy matrix
# Via seurat / cosine distance clustering
seuratSNN <- function(mat, resolution = 1, k.param = 10){ 
  set.seed(1)
  rownames(mat) <- make.unique(rownames(mat))
  obj <- FindNeighbors(mat, k.param = k.param, annoy.metric = "cosine")
  clusters <- FindClusters(object = obj$snn, resolution = resolution)
  return(as.character(clusters[,1]))
}

# Process new
archr_df <- readRDS(file = "../output/ArchR_main_metadata.rds")
archr_df$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(archr_df))
mSE <- readRDS("../output/mitoMutations_cov10_marrow.rds")
archr_cluster_vec <- archr_df$Clusters; names(archr_cluster_vec) <- as.character(archr_df$barcode) 
afp <- assays(mSE)[["allele_frequency"]]

# Call clones like in the mgatk/mtscatac paper
clonal_clusters <- seuratSNN(t(sqrt(afp)), resolution =  2.5, k.param = 10) 
clonal_clusters <- str_pad(as.character(as.numeric(factor(clonal_clusters))), 2, pad = "0")
table(clonal_clusters)
names(clonal_clusters) <- colnames(afp)

# Import protein data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,colnames(afp)]
adt_simple_norm <- data.matrix(t(t(adt_ss)/colSums(adt_ss))*1000)
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")

mat <- data.matrix(adtbm@assays$ADT@scale.data)

# Do a kruskal wallis test
protein_df <- data.frame(
  marker = rownames(mat)
)
protein_df$kruskal_pvalue <- sapply(1:dim(protein_df)[1], function(i){
  two_df <- data.frame(
    clocluster = clonal_clusters,
    protein =mat[i,] 
  )
  res.kw <- kruskal.test(protein ~ clocluster, data = two_df)
  (res.kw)[[3]]
})

protein_df$kruskal_perm_pvalue <- sapply(1:dim(protein_df)[1], function(i){
  set.seed(1)
  two_df <- data.frame(
    clocluster = sample(clonal_clusters),
    protein =mat[i,] 
  )
  res.kw <- kruskal.test(protein ~ clocluster, data = two_df)
  (res.kw)[[3]]
})
protein_df$kruskal_pvalue_adj <- p.adjust(protein_df$kruskal_pvalue)
rank_ordered <- protein_df %>% arrange((kruskal_pvalue)) 
rank_ordered$kruskal_perm_pvalue <- sort(rank_ordered$kruskal_perm_pvalue, decreasing = FALSE)
rank_ordered$rank <- 1:dim(rank_ordered)[1]

ggplot(rank_ordered, aes(x = rank, y = -log10(kruskal_pvalue_adj), color = kruskal_pvalue_adj < 0.0001)) + 
  geom_point(size = 0.3) + scale_color_manual(values = c("black", "firebrick"))  + 
  labs(x = "Rank sorted proteins", y = "-log10 Kruskal Wallis") + 
  pretty_plot(fontsize = 8) + L_border() + 
  theme(legend.position = "none")

library(Seurat)

mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
mat <- data.matrix(adtbm@assays$ADT@scale.data)

corPS <- cor(mdf$myeloidPS, t(mat), use = "pairwise.complete", method = "spearman")
corPS2 <- cor(mdf$myeloidPS2, t(mat), use = "pairwise.complete", method = "spearman")
cor_pvalue <- sapply(1:dim(mat)[1], function(i){
  cor.test(mdf$myeloidPS, mat[i,], use = "pairwise.complete", method = "spearman")$p.value
})

data.frame(ADT = rownames(mat),
           corPS = corPS[1,],
           corPS2 = corPS2[1,],
           pvalue = cor_pvalue) %>% arrange(desc(corPS)) -> cor_df

ggplot(cor_df,aes(x = corPS, y = -log10(pvalue), label = ADT)) +
  geom_text()

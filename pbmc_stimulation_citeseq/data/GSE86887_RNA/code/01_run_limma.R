library(ggplot2)
library(BuenColors)
library(limma)
library(dplyr)
library(data.table)
library(edgeR)
library(ggrastr)
"%ni%" <- Negate("%in%")

# Import gene expression
files <- list.files("../data/", full.names = TRUE)
counts <- sapply(files, function(x){
  fread(x)[[2]]
})
boo <- rowMeans(counts) > 10
counts <- counts[boo,]

# update names
genes <- fread(files[[1]])[[1]][boo]
rownames(counts) <- genes

# model with voom
design <- model.matrix(~ c(rep(1, 6), rep(0, 6)))
nf <-  calcNormFactors(counts)
y <- voom(counts = counts, design, plot = TRUE)

# Model
fit <- lmFit(y, design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
df <- topTable(fit, number = dim(counts)[1])
df$gene <- rownames(df)
df %>% filter(gene == "IFNG")
df %>% filter(logFC >2 & adj.P.Val < 0.00001) %>% pull(gene) %>% data.frame() -> activated_gene_df
write.table(activated_gene_df, file = "../output/CD3CD28_activated_gene_module.tsv", 
            sep = "", quote = FALSE, row.names = FALSE, col.names = FALSE)

vP <- ggplot(df, aes(x = -1*logFC, y = -1*log10(adj.P.Val))) + 
  geom_point_rast(size = 0.2, raster.dpi = 500) +
  pretty_plot(fontsize = 6) + labs(x = "logFC", y = "-log10 adj.P.Val") +
  L_border()
cowplot::ggsave2(vP, filename = "../output/differentialExpression_RNAseq_activated.pdf", 
                width = 2, height = 2)


library(pheatmap)
cor_mat <- cor(t(tfs), t(mat), use = "pairwise.complete")

pheatmap(cor_mat)
rv  <- rowVars(tfs)
threshold <- sort(rv, decreasing = TRUE)[150]
pheatmap(cor_mat[(rv > threshold) & (rowSums(is.na(cor_mat)) == 0),])

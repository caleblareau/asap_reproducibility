library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

call_mutations_mgatk <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10){
  
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  
  # Process mutation for one alternate letter
  process_letter <- function(letter){
    print(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    
    # Single cell functions
    getMutMatrix <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_fw  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_rev  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    # Bulk functions
    getBulk <- function(letter){
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
      return(vec)
    }
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    }
    
    update_missing_w_zero <- function(vec){
      ifelse(is.na(vec)  | is.nan(vec), 0, vec)
    }
    # Set up correlation per non-zero mutation based on the strands
    dt <- merge(data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo,])), 
                data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo,])), 
                by.x = c("i", "j"), by.y = c("i", "j"), 
                all = TRUE)[x.x >0 | x.y >0]
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    
    dt2 <- data.table(variant = variant_name[dt[[1]]],
                      cell_idx = dt[[2]], 
                      forward = dt[[3]],
                      reverse = dt[[4]])
    rm(dt)
    cor_dt <- dt2[, .(cor = cor(c(forward), c(reverse), method = "pearson", use = "pairwise.complete")), by = list(variant)]
    
    # Put in vector for convenience
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character(cor_dt$variant )
    
    # Compute the single-cell data
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i,dim(mat)[1]),
      j = c(summary(mat)$j,dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    
    # Compute bulk statistics
    mean = update_missing_w_zero(getBulk(letter))
    
    # Stablize variances by replacing low coverage cells with mean
    if(stabilize_variance){
      
      # Get indices of cell/variants where the coverage is low and pull the mean for that variant
      idx_mat <- which(data.matrix(cov[boo,] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[,1]]
      
      # Now, make sparse matrices for quick conversion
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = 1
      )
      
      means_mat <- sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      
      mmat2 <- mmat * ones + means_mat
      variance = rowVars(mmat2)
      rm(mmat2); rm(ones); rm(means_mat); rm(idx_mat); rm(idx_mat_mean)
      
    } else {
      variance = rowVars(mmat)
    }
    
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo,] >= 2) + (assays(SE)[[paste0(letter, "_counts_rev")]][boo,] >=2 )
    
    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide, 
      variant = variant_name,
      vmr = variance/(mean + 0.00000000001),
      mean = round(mean,7),
      variance = round(variance,7),
      n_cells_conf_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5 = Matrix::rowSums(mmat >= 0.05), 
      n_cells_over_10 = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20 = Matrix::rowSums(mmat >= 0.20),
      strand_correlation = cor_vec_val[variant_name],
      mean_coverage = Matrix::rowMeans(cov)[boo], 
      stringsAsFactors = FALSE, row.names = variant_name
    )
    se_new <- SummarizedExperiment(
      rowData = var_summary_df, 
      colData = colData(SE), 
      assays = list(allele_frequency = mmat, coverage = cov[boo,])
    )
    return(se_new)
  }
  
  return(SummarizedExperiment::rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
  
}
ctrl <- readRDS("LLL_CTRL_GRCh38-mtMask_atac_mgatk/final/LLL_CTRL_GRCh38-mtMask_atac_mgatk.rds") 
stim <- readRDS("LLL_STIM_GRCh38-mtMask_atac_mgatk/final/LLL_STIM_GRCh38-mtMask_atac_mgatk.rds")
table(stim$depth > 20)
table(ctrl$depth > 20)
stim_mgatk <- call_mutations_mgatk(stim[,stim$depth > 20])
ctrl_mgatk <- call_mutations_mgatk(ctrl[,ctrl$depth > 20])
both_mgatk <- call_mutations_mgatk(cbind(ctrl[,ctrl$depth > 20], stim[,stim$depth > 20]))

variants <- data.frame(rowData(both_mgatk)) %>%  dplyr::filter(n_cells_conf_detected >= 5 & strand_correlation > 0.65 & log10(vmr) > -2 ) %>% pull(variant)
ddf <- data.frame(
  variants = variants,
  ctrl_af = rowData(ctrl_mgatk)[variants,c("mean")],
  stim_af = rowData(stim_mgatk)[variants,c("mean")]
)
  
ggplot(ddf, aes(x = ctrl_af*100, y = stim_af*100, label = variants)) +
  geom_text() + scale_y_log10() + scale_x_log10() + pretty_plot() + L_border() +
  labs(x = "Control AF%", y = "Stim AF%")

cor(log10(ddf$ctrl_af), log10(ddf$stim_af))
data.frame(rowData(both_mgatk))[variants,] %>% arrange((mean))

mmdf <- merge(
  fread("LLL_CTRL_GRCh38-mtMask_atac_mgatk/final/LLL_CTRL_GRCh38-mtMask_atac_mgatk.depthTable.txt"),
  fread("LLL_CTRL_GRCh38-mtMask_rna_mgatk/final/LLL_CTRL_GRCh38-mtMask_rna_mgatk.depthTable.txt"),
  by = "V1"
)


ggplot(mmdf %>% dplyr::filter(V2.x < 250), aes(x = V2.x, y = V2.y)) + geom_point(alpha = 0.3) +
  pretty_plot() + L_border() + scale_y_log10() + scale_x_log10() + labs(x = "mean mtDNA coverage - ATAC", y = "mean mtDNA coverage - RNA")


posdf <- data.frame(
  pos = 1:16569,
  ctrl = rowMeans(assays(ctrl)[["coverage"]]),
  stim = rowMeans(assays(stim)[["coverage"]])
)
reshape2::melt(posdf , id.vars = "pos") %>%
  dplyr::filter(value > 20)%>% 
  ggplot(aes(x = pos, y = value, color = variable)) +
  geom_line() +
  pretty_plot() + L_border() +
  

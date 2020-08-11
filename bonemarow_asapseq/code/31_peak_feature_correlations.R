library(Signac)
library(Seurat)
library(ComplexHeatmap)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(data.table)
"%ni%" <- Negate("%in%")

# Get coordinates
gene.coords <- genes(EnsDb.Hsapiens.v86, filter = ~ gene_biotype == "protein_coding")
seqlevelsStyle(gene.coords) <- 'UCSC'
genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
extended_gene_cords <- Extend(x = gene.coords, upstream = 100000, downstream = 100000)
genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 2000, downstream = 0)

dat <- fread("../../multiple_datasets/data/marker_gene_mapping.tsv")

# Import the HQ cells
mdf <- readRDS("../output/ArchR_main_metadata.rds")
mdf$barcode <- gsub("ASAP_marrow_hg38#", "", rownames(mdf))

# Look at the ADT data
adt_ss <- readRDS("../output/adt_mat/ASAP_bonemarrow_matrix.rds")[,mdf$barcode]
adtbm <- CreateSeuratObject(counts = adt_ss, assay = "ADT")
adtbm <- NormalizeData(adtbm, assay = "ADT", normalization.method = "CLR")
adtbm <- ScaleData(adtbm, assay="ADT")
mat <- data.matrix(adtbm@assays$ADT@scale.data)

gs_mat <- readRDS("../../../asap_large_data_files/bonemarrow_data/output/signac_marrow_gene_scores.rds")

# Make a master data frame
master_df <- data.frame(
  dat, 
  adt.idx = match(dat$Marker_name,rownames(mat)),
  gs.idx = match(dat$Gene_symbol, rownames(gs_mat)),
  grange.idx = match(dat$Gene_symbol,extended_gene_cords$gene_name)
)
master_df_filt <- master_df[complete.cases(master_df),]

# Eventually won't need this but roll 
# get cells
cells <- intersect(colnames(gs_mat), colnames(mat))
mat <- mat[,cells]
gs_mat <- gs_mat[,cells]

# Import peaks x cells and filter for overlap with what's in the gene score / adt pairs
h5_raw <- Read10X_h5("../../../asap_large_data_files/bonemarrow_data/input/asap_marrow_hg38_raw_peak_bc_matrix.h5")
peak_10x_gr <- data.frame(str_split_fixed(rownames(h5_raw), "-|:", 3)) %>% setNames(c("chr", "start", "end")) %>% makeGRangesFromDataFrame()
ov <- findOverlaps(peak_10x_gr, extended_gene_cords[master_df_filt$grange.idx])
idx_keep <- unique(queryHits(ov))
peak_mat_filt <- h5_raw[idx_keep,cells]
peak_10x_filt_gr <- peak_10x_gr[idx_keep]
rm(h5_raw)

# Set up permutaiton
set.seed(1234)
perm_cells_idx <- sample(1:length(cells))


# Look for regulatory - feature correlations
lapply(1:dim(master_df_filt)[1], function(idx_loop){
  print(idx_loop)
  # Extract one marker / gene
  one <- master_df_filt[idx_loop,]
  adt_vec <- mat[one$adt.idx,]
  gs_vec <- gs_mat[one$gs.idx,]
  
  # Find only peaks in distal regulatory regions-- not those that are in the part of the gene score calculation
  ov_peaks_extended <- findOverlaps(peak_10x_filt_gr, extended_gene_cords[one$grange.idx])
  ov_peaks_gs <- findOverlaps(peak_10x_filt_gr, genebodyandpromoter.coords[one$grange.idx])
  possible <- unique(queryHits(ov_peaks_extended))
  extended_only_idx_peaks <- possible[possible %ni% queryHits(ov_peaks_gs)]
  
  chromatin <- data.matrix(peak_mat_filt[extended_only_idx_peaks,])
  
  if(dim(chromatin)[1] > 0){
    # Now loop over applicable chromatin peaks
    lapply(1:length(rownames(chromatin)), function(peak_idx){
      
      # Do a bunch of correlation tests
      real_adt <- cor.test(chromatin[peak_idx,], adt_vec)
      permuted_adt <- cor.test(chromatin[peak_idx,], adt_vec[perm_cells_idx])
      real_gs <- cor.test(chromatin[peak_idx,], gs_vec)
      permuted_gs <- cor.test(chromatin[peak_idx,], gs_vec[perm_cells_idx])
      
      data.frame(one[,c(1,2)], peak = rownames(chromatin)[peak_idx], 
                 pvalue_adt = real_adt$p.value, cor_adt = real_adt$estimate,
                 pvalue_adt_perm = permuted_adt$p.value, cor_adt_perm = permuted_adt$estimate,
                 pvalue_gs = real_gs$p.value, cor_gs = real_gs$estimate,
                 pvalue_gs_perm = permuted_gs$p.value, cor_gs_perm = permuted_gs$estimate
      )
    }) %>% rbindlist() %>% data.frame() -> all_peaks_df
    #all_peaks_df$upstream <- upstream
    all_peaks_df
  } else {
    NULL
  }
}) %>% rbindlist() %>% data.frame() -> assoc_df

ggplot(assoc_df, aes(x = cor_adt, y = cor_gs)) + 
  geom_point() + labs(x = "Correlation ADT", y = "Correlation Gene Score") +
  pretty_plot() + L_border() + 
  ggtitle("Each dot is 1 peak - feature correlation")

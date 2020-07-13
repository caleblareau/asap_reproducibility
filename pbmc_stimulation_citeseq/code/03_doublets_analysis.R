library(data.table)
library(dplyr)
library(Matrix)
library(Seurat)

# Import scRNA-seq
import_scRNAseq <- function(dir_base){
  
  data.dir <- paste0("../data/rnaseq/", dir_base)
  raw <- Read10X(data.dir = data.dir)
  colnames(raw) <- paste0(substr(colnames(raw), 1, 16), "-1")
  
  # import scrublet results
  singlets <- fread(paste0("../data/rnaseq/scrublet_out/", dir_base, ".scrub.tsv")) %>%
    data.frame() %>% dplyr::filter(score < 100) %>% pull(barcode) # the original called threshold seemed too conservative; this is a better estimate for these libraries
  
  # Filter for singlet cells and non-mito genes
  raw <- raw[!grepl("^MT", rownames(raw)), singlets]
  raw
}

# Import ADT
import_kite_counts <- function(library){
  mtx <- fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0("../data/adt/",library,"_featurecounts/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  return(t(matx))
}

# Import matrices
ctrl_scRNA <- import_scRNAseq("ctrl")
stim_scRNA <- import_scRNAseq("stim")
ctrl_ADT <- import_kite_counts("ctrl")[,gsub("-1", "", colnames(ctrl_scRNA))]
stim_ADT <- import_kite_counts("stim")[,gsub("-1", "", colnames(stim_scRNA))]

booctrl <- !(colSums(ctrl_ADT) > 25000)
boostim <- !(colSums(stim_ADT) > 25000)

ctrl.cite <- CreateSeuratObject(counts = ctrl_ADT[,booctrl], assay = "ADT")
ctrl.cite <- NormalizeData(ctrl.cite, assay = "ADT", normalization.method = "CLR")
ctrl.cite <- ScaleData(ctrl.cite, assay = "ADT")

stim.cite <- CreateSeuratObject(counts = stim_ADT[,boostim], assay = "ADT")
stim.cite <- NormalizeData(stim.cite, assay = "ADT", normalization.method = "CLR")
stim.cite <- ScaleData(stim.cite, assay = "ADT")

rdf <- rbind(
  data.frame(t(ctrl.cite@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "CITE", what = "Control"),
  data.frame(t(stim.cite@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "CITE", what = "Stim"))

ggplot(rdf, aes(x = CD11c, y = CD4.1)) +
  facet_wrap(~what) + geom_point(alpha = 0.5)


process_rdf_props <- function(rdf, what, assay = "CITE"){
  rdf %>%
    mutate(pCD20 = CD20 > 2, pCD4 = CD4.1 > 0, pCD8 = CD8a > 1, pCD11c = CD11c > 1.8) -> ppdf
  
  denom <- dim(ppdf)[1]/100
  data.frame(
    what, assay, total_cells = denom*100,
    B_CD4 = sum(ppdf$pCD20 & ppdf$pCD4)/denom,
    B_CD8 = sum(ppdf$pCD20 & ppdf$pCD8)/denom,
    Mono_CD4 = sum(ppdf$pCD11c & ppdf$pCD4)/denom,
    Mono_CD8 = sum(ppdf$pCD11c & ppdf$pCD8)/denom,
    B_Mono = sum(ppdf$pCD20 & ppdf$pCD11c)/denom,
    CD8_CD4 = sum(ppdf$pCD8 & ppdf$pCD4)/denom
  )
}

rdf %>% filter(what == "Stim") %>% process_rdf_props(what = "Stim", assay = "CITE")
rdf %>% filter(what == "Control") %>% process_rdf_props(what = "Control", assay = "CITE")

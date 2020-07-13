library(data.table)
library(dplyr)
library(Seurat)
library(Matrix)
library(BuenColors)
library(harmony)

# Process kite counts
import_kite_counts <- function(library, bc_file, tech, bio){
  
  # Import the goodies
  bcs <- substr(fread(bc_file, header = FALSE)[[1]], 1, 16)
  mtx <- fread(paste0(library,"featurecounts.mtx.gz"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  
  # Label features
  rownames(matx) <- fread(paste0(library,"/featurecounts.barcodes.txt.gz"), header = FALSE)[[1]]
  colnames(matx) <- paste0(fread(paste0(library,"/featurecounts.genes.txt.gz"), header = FALSE)[[1]])
  maty <- t(matx)[,rownames(matx) %in% bcs]
  maty <- maty[,Matrix::rowSums(maty) >= 500 & Matrix::rowSums(maty) < 100000]
  colnames(maty) <- paste0(library, "_", colnames(maty))
  
  
  # Now create a seurat object
  raw <- CreateSeuratObject(counts = maty,  min.cells = 3, min.features = 10, assay = "ADT"); raw$tech <- tech; raw$bio <- bio
  raw <- NormalizeData(raw, assay = "ADT", normalization.method = "CLR")
  raw <- ScaleData(raw, assay = "ADT")
  
  return(raw)
}

# Import counts 
asap_control <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_ctrl_ADT/",  "../../pbmc_stimulation_asapseq/data/cellranger/control_barcodes.tsv", "ASAP", "noStim")
asap_stim <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_stim_ADT/",  "../../pbmc_stimulation_asapseq/data/cellranger/stim_barcodes.tsv","ASAP", "Stim")
cite_control <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/ctrl_featurecounts/",  "../../pbmc_stimulation_citeseq/data/rnaseq/ctrl/barcodes.tsv.gz","CITE", "noStim")
cite_stim <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/stim_featurecounts/",  "../../pbmc_stimulation_citeseq/data/rnaseq/stim/barcodes.tsv.gz","CITE", "Stim")

rdf <- rbind(
  data.frame(t(asap_control@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "ASAP", what = "Control"),
  data.frame(t(asap_stim@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "ASAP", what = "Stim"),
  data.frame(t(cite_control@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "CITE", what = "Control"),
  data.frame(t(cite_stim@assays$ADT@scale.data[c("CD19", "CD20", "CD11c", "CD14", "CD16", "CD8a",  "CD4-1", "CD3-1"),]),
             assay = "CITE", what = "Stim"))

ggplot(rdf %>% filter(CD20 < 8), aes(x = CD11c, y = CD20)) +
  facet_grid(assay~what) + geom_point(alpha = 0.5) +
  pretty_plot() +
  geom_vline(xintercept = 1.8, color = "firebrick") +
  geom_hline(yintercept = 2, color = "firebrick")



process_rdf_props <- function(rdf_in, what1, assay1){
  rdf_in  %>% filter(what == what1 & assay == assay1) %>%
    mutate(pCD20 = CD20 > 2, pCD4 = CD4.1 > 0, pCD8 = CD8a > 0.5, pCD11c = CD11c > 1.8) -> ppdf
  
  denom <- dim(ppdf)[1]/100
  data.frame(
    what1, assay1, total_cells = denom*100,
    B_CD4 = round(sum(ppdf$pCD20 & ppdf$pCD4)/denom,2),
    B_CD8 = round(sum(ppdf$pCD20 & ppdf$pCD8)/denom,2),
    Mono_CD4 = round(sum(ppdf$pCD11c & ppdf$pCD4)/denom,2),
    Mono_CD8 = round(sum(ppdf$pCD11c & ppdf$pCD8)/denom,2),
    B_Mono = round(sum(ppdf$pCD20 & ppdf$pCD11c)/denom,2),
    CD8_CD4 = round(sum(ppdf$pCD8 & ppdf$pCD4)/denom,2)
  )
}

rbind(rdf  %>% process_rdf_props(what1 = "Stim", assay1 = "ASAP"),
      rdf %>% process_rdf_props(what1 = "Control", assay1 = "ASAP"),
      rdf  %>% process_rdf_props(what1 = "Stim", assay1 = "CITE"),
      rdf %>% process_rdf_props(what1 = "Control", assay1 = "CITE")
)
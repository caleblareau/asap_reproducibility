library(data.table)

import_kite_counts <- function(path,library, gz = FALSE){
  if(gz){
    xx <- ".gz"
  } else {
    xx <- ""
  }
  mtx <- fread(paste0(path, "featurecounts",library,".mtx",xx), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".barcodes.txt",xx), header = FALSE)[[1]], "-1")
  colnames(matx) <- paste0(fread(paste0(path, "featurecounts",library,".genes.txt",xx), header = FALSE)[[1]])
  pct_human <- (matx[,2])/rowSums(matx)
  return(data.frame(barcode = rownames(matx), data.matrix(matx)))
}

path <- "../../../../asap_paper/species_mixing_raw_organize/"
sc1 <- merge(fread(paste0(path, "singlecell1.csv.gz")), import_kite_counts(path, "1"), by = "barcode")
sc2 <- merge(fread(paste0(path, "singlecell2.csv.gz")), import_kite_counts(path, "2"), by = "barcode")
write.table(sc1, file = "../output/SpeciesMix1_QC.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(sc2, file = "../output/SpeciesMix2_QC.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#----- Marrow ----
adt_bm <- import_kite_counts("../../bonemarow_asapseq/data/adt/", "", TRUE)
hto_bm <- import_kite_counts("../../bonemarow_asapseq/data/hto/", "", FALSE)
write.table(adt_bm, file = "../output/Human_BoneMarrow_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(hto_bm, file = "../output/Human_BoneMarrow_HTO.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#----- PBMC stim ----
a_c <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_ctrl_ADT/", "", TRUE)
a_s <- import_kite_counts("../../pbmc_stimulation_asapseq/data/adt/ASAP_stim_ADT/", "", TRUE)
c_c <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/ctrl_featurecounts/", "", TRUE)
c_s <- import_kite_counts("../../pbmc_stimulation_citeseq/data/adt/stim_featurecounts/", "", TRUE)

write.table(a_c, file = "../output/CD28_CD3_control_ASAP_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(a_s, file = "../output/CD28_CD3_stim_ASAP_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(c_c, file = "../output/CD28_CD3_control_CITE_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(c_s, file = "../output/CD28_CD3_stim_CITE_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#----- Broad PBMCs ----
broad_B <- import_kite_counts("../../broad_pbmc_asapseq/data/tag_data_analyses/tag_data_LL/ADT3_ASAP_B/", "", FALSE)
broad_C <- import_kite_counts("../../broad_pbmc_asapseq/data/tag_data_analyses/tag_data_LL/HTO1_ASAP_C/", "", FALSE)
write.table(broad_B[,c(1:10)], file = "../output/Broad_LibB_ADT.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(broad_C[,c(1, 11:14)], file = "../output/Broad_LibC_HTO.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


#----- NYGC PBMCs ----
intraA <- import_kite_counts("../../nygc_pbmcs/data/intra/", "A", TRUE)
intraB <- import_kite_counts("../../nygc_pbmcs/data/intra/", "B", TRUE)
LLLA <- import_kite_counts("../../nygc_pbmcs/data/AvsB/TSA_LLL/featurecounts/", "", FALSE)
LLLB <- import_kite_counts("../../nygc_pbmcs/data/AvsB/TSB_LLL/featurecounts/", "", FALSE)
OA <- import_kite_counts("../../nygc_pbmcs/data/AvsB/TSA_OMNI/featurecounts/", "", FALSE)
OB <- import_kite_counts("../../nygc_pbmcs/data/AvsB/TSB_OMNI/featurecounts/", "", FALSE)

write.table(LLLA, file = "../output/AvsB_TSA_LLL.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(LLLB, file = "../output/AvsB_TSB_LLL.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(OA, file = "../output/AvsB_TSA_OMNI.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(OB, file = "../output/AvsB_TSB_OMNI.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(intraA, file = "../output/Intra_TSA.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(intraB, file = "../output/Intra_TSB.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)





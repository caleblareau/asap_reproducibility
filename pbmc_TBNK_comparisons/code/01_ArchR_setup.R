library(data.table)
library(ArchR)
addArchRGenome("hg38")


arrow_omni <- createArrowFiles(
  inputFiles = "../../../asap_large_data_files/nygc_pbmc/input/AvsB_ATAC_OMNI_fragments.tsv.gz",
  sampleNames = "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_OMNI",
  filterTSS = 4, 
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

arrow_LLL <- createArrowFiles(
  inputFiles = "../../../asap_large_data_files/nygc_pbmc/input/AvsB_ATAC_LLL_fragments.tsv.gz",
  sampleNames = "../../../asap_large_data_files/nygc_pbmc/output/archr/AvsB_LLL",
  filterTSS = 4, 
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

arrow_intra <- createArrowFiles(
  inputFiles = "../../../asap_large_data_files/nygc_pbmc/input/Intra_ATAC_fragments.tsv.gz",
  sampleNames = "PBMCs_Intra",
  filterTSS = 4, 
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)


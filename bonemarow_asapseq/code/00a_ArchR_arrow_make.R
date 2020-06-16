library(ArchR)

addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = "../../../asap_large_data_files/bonemarrow_data/input",
  sampleNames = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38",
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

proj <- ArchRProject(
  ArrowFiles = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38.arrow", 
  outputDirectory = "../../../asap_large_data_files/bonemarrow_data/output/archr_marrow/ASAP_marrow_hg38",
  copyArrows = FALSE
)

df <- data.frame(barcode = substr(rownames(proj@cellColData), 18, 35))
write.table(df, file = "../data/barcodes/step1_barcodes_from_ArchR.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
library(data.table)
library(Matrix)
library(ggrastr)

import_kite_counts <- function(path, library){
  mtx <- fread(paste0(path, "featurecounts",library,".mtx"), header = FALSE)
  dim <- mtx[1,]
  mtx <- mtx[-1,]
  matx <- sparseMatrix(i = mtx[[1]], j = mtx[[2]], x = mtx[[3]])
  rownames(matx) <- paste0(fread(paste0(path,"featurecounts",library,".barcodes.txt"), header = FALSE)[[1]], "-1")
  colnames(matx) <- paste0(fread(paste0(path,"featurecounts",library,".genes.txt"), header = FALSE)[[1]])
  pct_human <- (matx[,2])/Matrix::rowSums(matx)
  return(data.frame(barcode = rownames(matx), data.matrix(matx),
                    total_cite = Matrix::rowSums(matx),
                    pct_human_cite = pct_human*100))
}

assign_organism <- function(df){
  df$pct_human_atac <- (df$peak_region_fragments_GRCh38)/(df$peak_region_fragments_GRCh38 + df$peak_region_fragments_mm10) * 100
  df$pct_human_cite <- df$pct_human_cite
  
  pha <- df$pct_human_atac
  phc <- df$pct_human_cite
  dha <- df$peak_region_fragments_GRCh38 > 100
  dhc <- df$hCD29 > 100
  dma <- df$peak_region_fragments_mm10 > 100
  dmc <- df$mCD29 > 50
  
  df$assignment_atac <- case_when(
    dha & (pha > 95) ~ "human",
    dma & (pha < 5) ~ "mouse",
    dha & dma ~ "mixed",
    TRUE ~ "none"
  )
  
  df$assignment_cite <- case_when(
    dhc & (phc > 90) ~ "human",
    dmc & (phc < 10) ~ "mouse",
    dhc & dmc ~ "mixed",
    TRUE ~ "none"
  )
  
  df
}

one <- import_kite_counts("../data/adt/", "1")
two <- import_kite_counts("../data/adt/", "2")
sc_one <- merge(one, fread("../data/singlecell1.csv.gz"), by = "barcode") %>% assign_organism
sc_two <- merge(two, fread("../data/singlecell2.csv.gz"), by = "barcode") %>% assign_organism

p1 <- ggplot(sc_one , aes(x = peak_region_fragments_GRCh38/1000, y = peak_region_fragments_mm10/1000, color = assignment_atac)) +
  geom_point_rast(size = 2)  + pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "purple2","firebrick", "lightgrey")) + 
  theme(legend.position = "none") +
  labs(y = "mm10 fragments in peaks", x = "GRCh38 fragments in peaks", color = "ATAC assign")

p2 <- ggplot(sc_one , aes(y = mCD29, x = hCD29, color = assignment_atac)) +
  geom_point_rast(size = 2) + pretty_plot(fontsize = 8) + L_border() + 
  scale_x_continuous(limits = c(0, 2500)) + scale_y_continuous(limits = c(0, 2500)) +
  scale_color_manual(values = c("dodgerblue3", "purple2","firebrick", "lightgrey")) + 
  theme(legend.position = "none") +
  labs(y = "mCD29 ADT UBIs", x = "hCD29 ADT UBIs", color = "ATAC assign")

cowplot::ggsave2(cowplot::plot_grid(p1, p2, nrow =1), width = 3.6, height = 1.8, file = "../plots/mixing_plots_preSPRI.pdf")

p1t <- ggplot(sc_two , aes(x = peak_region_fragments_GRCh38/1000, y = peak_region_fragments_mm10/1000, color = assignment_atac)) +
  geom_point_rast(size = 2)  + pretty_plot(fontsize = 8) + L_border() + 
  scale_color_manual(values = c("dodgerblue3", "purple2","firebrick", "lightgrey")) + 
  theme(legend.position = "none") +
  labs(y = "mm10 fragments in peaks", x = "GRCh38 fragments in peaks", color = "ATAC assign")

p2t <- ggplot(sc_two , aes(y = mCD29, x = hCD29, color = assignment_atac)) +
  geom_point_rast(size = 2) + pretty_plot(fontsize = 8) + L_border() + 
  scale_x_continuous(limits = c(0, 2500)) + scale_y_continuous(limits = c(0, 2500)) +
  scale_color_manual(values = c("dodgerblue3", "purple2","firebrick", "lightgrey")) + 
  theme(legend.position = "none") +
  labs(y = "mCD29 ADT UBIs", x = "hCD29 ADT UBIs", color = "ATAC assign")

cowplot::ggsave2(cowplot::plot_grid(p1t, p2t, nrow =1), width = 3.6, height = 1.8, file = "../plots/mixing_plots_postSPRI.pdf")





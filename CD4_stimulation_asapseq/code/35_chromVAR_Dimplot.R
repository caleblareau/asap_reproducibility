options(stringsAsFactors=F)

#******************************************
#### args        = commandArgs(trailingOnly = T)
#### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP="Perturb_CD4_stim"

#******************************************
suppressMessages(library("patchwork"))
suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrastr"))

dir.create("plots/chromVAR_motif_Dimplot", showWarnings = FALSE)

#******************************************
##### data 
#******************************************

UMAP    =fread("../output/ArchR/DimPlot/ArchR_UMAP_res.txt") %>% as.data.frame(.)

chromVAR_res = fread("../output/ArchR/chromVAR/motif_zscore_ChromVAR.txt") %>% as.data.frame(.) %>% column_to_rownames("motif") %>% t() %>% as.data.frame(.) %>%
                rownames_to_column("UMI")

chromVAR_dimPlot = full_join(UMAP,chromVAR_res,by="UMI")
colnames(chromVAR_dimPlot)[-c(1:3)] = colnames(chromVAR_dimPlot)[-c(1:3)] %>% strsplit(.,"_") %>%
                                       sapply(.,function(x){x[1]})

#### plot figures using ggplot2
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
col_list   = colorRampPalette(blueYellow)

###########
motif_tmp = "BACH1"
min_score = -5
max_score = 5
tmp = chromVAR_dimPlot[,motif_tmp]
tmp[tmp<min_score] = min_score
tmp[tmp>max_score] = max_score

chromVAR_dimPlot_tmp = chromVAR_dimPlot
chromVAR_dimPlot_tmp[,motif_tmp] = tmp

p1 = ggplot(chromVAR_dimPlot_tmp, aes_string(x="UMAP_1",y="UMAP_2",color=motif_tmp)) + 
      geom_point_rast(size = .5) + 
      theme_void() + 
      scale_color_gradientn(colours=col_list(250))+
      labs(x="UMAP_1",y="UMAP_2",title = paste0(SAMPLE_TMP," ",motif_tmp," chromVAR deviation zscore"))

paste0("../plots/chromVAR_motif_Dimplot/chromVAR_deviation_zscore_",motif_tmp,"_Dimplot.pdf") %>% pdf(.,h=5,w=5)
 plot(p1)
dev.off()

###########
motif_tmp = "ETS1"
min_score = -2
max_score = 2
tmp = chromVAR_dimPlot[,motif_tmp]
tmp[tmp<min_score] = min_score
tmp[tmp>max_score] = max_score

chromVAR_dimPlot_tmp = chromVAR_dimPlot
chromVAR_dimPlot_tmp[,motif_tmp] = tmp

p1 = ggplot(chromVAR_dimPlot_tmp, aes_string(x="UMAP_1",y="UMAP_2",color=motif_tmp)) + 
      geom_point_rast(size = .5) + 
      theme_void() + 
      scale_color_gradientn(colours=col_list(250))+
      labs(x="UMAP_1",y="UMAP_2",title = paste0(SAMPLE_TMP," ",motif_tmp," chromVAR deviation zscore"))

paste0("../plots/chromVAR_motif_Dimplot/chromVAR_deviation_zscore_",motif_tmp,"_Dimplot.pdf") %>% pdf(.,h=5,w=5)
 plot(p1)
dev.off()

###########
motif_tmp = "NFKB1"
min_score = -4
max_score = 4
tmp = chromVAR_dimPlot[,motif_tmp]
tmp[tmp<min_score] = min_score
tmp[tmp>max_score] = max_score

chromVAR_dimPlot_tmp = chromVAR_dimPlot
chromVAR_dimPlot_tmp[,motif_tmp] = tmp

p1 = ggplot(chromVAR_dimPlot_tmp, aes_string(x="UMAP_1",y="UMAP_2",color=motif_tmp)) + 
      geom_point_rast(size = .5) + 
      theme_void() + 
      scale_color_gradientn(colours=col_list(250))+
      labs(x="UMAP_1",y="UMAP_2",title = paste0(SAMPLE_TMP," ",motif_tmp," chromVAR deviation zscore"))

paste0("../plots/chromVAR_motif_Dimplot/chromVAR_deviation_zscore_",motif_tmp,"_Dimplot.pdf") %>% pdf(.,h=5,w=5)
 plot(p1)
dev.off()

###########
motif_tmp = "NFATC2"
min_score = -4
max_score = 4
tmp = chromVAR_dimPlot[,motif_tmp]
tmp[tmp<min_score] = min_score
tmp[tmp>max_score] = max_score

chromVAR_dimPlot_tmp = chromVAR_dimPlot
chromVAR_dimPlot_tmp[,motif_tmp] = tmp

p1 = ggplot(chromVAR_dimPlot_tmp, aes_string(x="UMAP_1",y="UMAP_2",color=motif_tmp)) + 
      geom_point_rast(size = .5) + 
      theme_void() + 
      scale_color_gradientn(colours=col_list(250))+
      labs(x="UMAP_1",y="UMAP_2",title = paste0(SAMPLE_TMP," ",motif_tmp," chromVAR deviation zscore"))

################ Fig E ################################################################################
paste0("../plots/chromVAR_motif_Dimplot/chromVAR_deviation_zscore_",motif_tmp,"_Dimplot.pdf") %>% pdf(.,h=5,w=5)
 plot(p1)
dev.off()

options(stringsAsFactors=F)

#******************************************
#### args        = commandArgs(trailingOnly = T)
#### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP="Perturb_CD4_stim"

#******************************************
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(pheatmap))

#################
sgRNA_list     = read.table("../data/sgRNA_list.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")
HTO_res = paste0("../output/Signac/after_filter_Signac/HTO_res_filtered.txt") %>% fread(.) %>% as.data.frame(.)
ADT_CLR = paste0("../output/Signac/after_filter_Signac/ADT_filtered.txt") %>% read.table(.,header=T,row.names=1,stringsAsFactors=F)
ADT_CLR_2 = ADT_CLR %>% rownames_to_column("CellBarcode") %>% full_join(HTO_res,.,by="CellBarcode")
colnames(ADT_CLR_2)[-c(1:3)] = colnames(ADT_CLR_2)[-c(1:3)] %>% paste0("ADT_",.)

ADT_CLR_3 = ADT_CLR_2 %>% group_by(HTO) %>% 
             summarise_each(funs(mean), starts_with("ADT_"))  %>% as.data.frame() %>% column_to_rownames("HTO")
colnames(ADT_CLR_3) = gsub("ADT_","",colnames(ADT_CLR_3))

###################

target_ADT = c("CD3","CD4","CD8","TCRAB","CD28","CD45","CD45RA","CD45RO","CD197.CCR7",
               "CD62L","CD95.Fas","CD27","CD31","CD38","CD127.IL.7R","CD25",
               "CD69","CD279.PD.1","OX40","CD137.4.1BB","CD357.GITR","CD366.Tim.3",
               "CD39","CD278.ICOS","CD7","TIGIT.VSTM3","CD73","Mouse.IgG1")
ADT_CLR_4 = ADT_CLR_3[rev(sgRNA_list $HTO),target_ADT]


############pheatmap##################

col=c("#6594C4", "#77A6CE", "#8BB9D7", "#9EC7DF", "#B2D4E7", "#C6E2EE", "#DAEEF5",
      "#E5F5EC", "#EDF8DE", "#F4FBD0", "#FDFDC1", "#FEF8B5", "#FEF0A9", "#FEE89D",
      "#FEE091", "#FDCE83", "#FDB976", "#FCA468", "#FC8F5A", "#F4794E")

################ Fig C ################################################################################
pheatmap(ADT_CLR_4, cluster_rows=FALSE, cluster_cols = FALSE, cellwidth = 20, cellheight = 20, color=col, 
         filename=paste0("../plots/ADT_CLR_mean_sgRNA_heatmap.pdf"))


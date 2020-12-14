options(stringsAsFactors=F)

#******************************************
### args        = commandArgs(trailingOnly = T)
### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP="Perturb_CD4_stim"

#******************************************
suppressMessages(library("patchwork"))
suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrastr"))
suppressMessages(library("BuenColors"))
suppressMessages(library("Seurat"))

#dir_list = c("output/ADT/ALLsgRNA_Dimplot","output/ADT/each_sgRNA_Dimplot",
#             "plots/ADT_ALLsgRNA_Dimplot","plots/ADT_each_sgRNA_Dimplot")
#for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}

#******************************************
##### data 
#******************************************

UMAP    = fread("../output/ArchR/DimPlot/ArchR_UMAP_res.txt") %>% as.data.frame(.)
UMAP$CellBarcode = UMAP$UMI %>% strsplit(.,"#") %>% sapply(.,function(x){x[2]})
ADT_CLR = read.table("../output/Signac/after_filter_Signac/ADT_filtered.txt",header=T,row.names=NULL,stringsAsFactors=F)
HTO_res = fread("../output/Signac/after_filter_Signac/HTO_res_filtered.txt") %>% as.data.frame(.)
HTO_sgRNAlist  = c("sgNTC", "sgCD4","sgCD3E","sgCD3ECD4","sgZAP70","sgNFKB2")

ADT_CLR_Dimplot = full_join(UMAP[,c("CellBarcode","UMAP_1","UMAP_2")],HTO_res,by="CellBarcode") %>% 
                   full_join(.,ADT_CLR,by="CellBarcode")

###### percentile set to 05-95 ######
percntile_adjust_up   = function(DATA,perc=0.90){thresh=quantile(DATA,perc);DATA[DATA>thresh]=thresh;DATA}
percntile_adjust_down = function(DATA,perc=0.05){thresh=quantile(DATA,perc);DATA[DATA<thresh]=thresh;DATA}

ADT_CLR_Dimplot_2 = ADT_CLR_Dimplot
for(TARGET_tmp in colnames(ADT_CLR_Dimplot)[6:ncol(ADT_CLR_Dimplot)]){
  ADT_CLR_Dimplot_2[,TARGET_tmp] =  ADT_CLR_Dimplot[,TARGET_tmp] %>% 
                                       percntile_adjust_up(.) %>%
                                       percntile_adjust_down(.)
 }

#******************************************
##### drwaing ADT Dimplot about all ADT data 
#******************************************
for(TARGET_tmp in colnames(ADT_CLR_Dimplot)[6:ncol(ADT_CLR_Dimplot)]){
  ADT_CLR_Dimplot_tmp = ADT_CLR_Dimplot_2 %>% dplyr::select("UMAP_1","UMAP_2","HTO_sgRNA",TARGET_tmp)
  colnames(ADT_CLR_Dimplot_tmp)[4] = "TARGET"
  ADT_CLR_Dimplot_tmp = ADT_CLR_Dimplot_tmp %>% arrange(TARGET)

  p1 = ggplot(ADT_CLR_Dimplot_tmp, aes(x=UMAP_1,y=UMAP_2,color=TARGET)) + 
        geom_point_rast(size = .5) + 
        theme_void() + 
        scale_color_gradientn(colours=jdb_palette("solar_basic"))+
        labs(x="UMAP_1",y="UMAP_2",title = paste0(SAMPLE_TMP," ",TARGET_tmp," ADT_CLR score"))

  paste0("../output/ADT/ALLsgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot.pdf") %>% pdf(.,h=5,w=5)
   plot(p1)
  dev.off()

  #####
  dim_Plot = list()
  for(tmp_num in 1:6){
      dim_Plot[[tmp_num]] = ggplot(ADT_CLR_Dimplot_tmp, aes(x=UMAP_1,y=UMAP_2))+
              geom_point(size = .5,col="gray85") +
              geom_point(data = ADT_CLR_Dimplot_tmp %>% dplyr::filter(HTO_sgRNA==HTO_sgRNAlist[tmp_num]),aes(x=UMAP_1,y=UMAP_2,col=TARGET),size = 1) +
              scale_color_gradientn(colours=jdb_palette("solar_basic"),limits=c(0,max(ADT_CLR_Dimplot_tmp$TARGET)))+
              theme_void() +
              labs(x="UMAP_1",y="UMAP_2")
      dim_Plot[[tmp_num]] = AugmentPlot(dim_Plot[[tmp_num]], h=5, w=5, dpi=100)
  }

   tmp_num=7
   dim_Plot[[tmp_num]] = ggplot(ADT_CLR_Dimplot_tmp, aes(x=UMAP_1,y=UMAP_2,col=TARGET))+
              geom_point(size = 1) +
              scale_color_gradientn(colours=jdb_palette("solar_basic"),limits=c(0,max(ADT_CLR_Dimplot_tmp$TARGET)))+
              theme_void() +
              labs(x="UMAP_1",y="UMAP_2")
   dim_Plot[[tmp_num]] = AugmentPlot(dim_Plot[[tmp_num]], h=5, w=5, dpi=100) + plot_annotation(title = "ALL")

  p2 = (dim_Plot[[1]]|dim_Plot[[2]]|dim_Plot[[3]]|dim_Plot[[4]]|dim_Plot[[5]]|dim_Plot[[6]]|dim_Plot[[7]]) + 
        plot_annotation(title = paste0(SAMPLE_TMP," ",TARGET_tmp," ADT_CLR score",paste0(HTO_sgRNAlist,collapse="_")))

  paste0("../output/ADT/each_sgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot_each_sgRNA.pdf") %>% pdf(.,h=5,w=35)
   plot(p2)
  dev.off()
}

#******************************************
##### copying data needed
#******************************************

############ SupFig C
for(TARGET_tmp in c("CD3","CD4","CD25")){
  file.copy(paste0("../output/ADT/each_sgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot_each_sgRNA.pdf"),
            paste0("../plots/ADT_each_sgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot_each_sgRNA.pdf") )
}

############ SupFig D
for(TARGET_tmp in c("CD197.CCR7","CD95.Fas","CD69","OX40",
                    "CD137.4.1BB","CD357.GITR","CD279.PD.1","CD366.Tim.3")){
  file.copy(paste0("../output/ADT/ALLsgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot.pdf"),
            paste0("../plots/ADT_ALLsgRNA_Dimplot/ADT_CLR_",TARGET_tmp,"_Dimplot.pdf") )
}



options(stringsAsFactors=F)

#******************************************
### args        = commandArgs(trailingOnly = T)
### SAMPLE_TMP  = as.character( args[1] )

SAMPLE_TMP  = "Perturb_CD4_stim"

##################
suppressMessages(library("ArchR"))
suppressMessages(library("patchwork"))
suppressMessages(library("tidyverse"))
suppressMessages(library("gridExtra"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrepel"))

set.seed(111)

#dir_list = c("output/ArchR/chromVAR/volcanoPlot", "plots/chromVAR_volcanoPlot","output/ArchR/chromVAR/scatterPlot_logFDR","output/ArchR/chromVAR/scatterPlot_median_diff",
#             "plots/chromVAR_scatterPlot_logFDR","plots/chromVAR_scatterPlot_median_diff")
#for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}

#****************************************#
# Loading information
#****************************************#
load("../../../asap_large_data_files/CD4_stimulation_asapseq/output/ArchR.Rdata")

#****************************************#
# summarizing WMW analysis
#****************************************#

med_res  = fread("../output/ArchR/chromVAR/WMW_res/deviation_zScore_median_data.txt") %>% as.data.frame(.)

for(sgRNA_TMP in sgRNA_list$HTO[3:10]){
  df_tmp        = paste0("../output/ArchR/chromVAR/WMW_res/",sgRNA_TMP,"_WMW_summary.txt") %>% 
                   fread(.) %>% as.data.frame(.)
  colnames(df_tmp) = c("sgRNA_TF","WMW_pValue")
  df_tmp$sgRNA  = strsplit(df_tmp$sgRNA_TF,"__") %>% sapply(.,function(x){x[1]})
  df_tmp$TF     = strsplit(df_tmp$sgRNA_TF,"__") %>% sapply(.,function(x){x[2]})
  df_tmp$WMW_qValue = p.adjust(df_tmp$WMW_pValue,method="BH")
  df_tmp$WMW_qValue[df_tmp$WMW_qValue==0] = 1/10^300         #################################################
  df_tmp$WMW_logFDR = -log10(df_tmp$WMW_qValue)

  med_res_tmp   = med_res[,c("TF","sgNTC",sgRNA_TMP)]
  med_res_tmp_2 = data.frame(TF=med_res_tmp$TF,diff_med=med_res_tmp[,3]-med_res_tmp[,2])

  df_tmp_2 = full_join(df_tmp,med_res_tmp_2,by="TF")
  df_tmp_2$Diff   = "non_significant"
  df_tmp_2$Diff[(df_tmp_2$WMW_qValue<0.05)&(abs(df_tmp_2$diff_med)>0.25)] = "significant"

  if(sgRNA_TMP==sgRNA_list$HTO[3]){WMW_sum=c()}
  WMW_sum = rbind(WMW_sum,df_tmp_2)
}

WMW_sum$sign_WMW_logFDR = (sapply(WMW_sum$diff_med,function(x){if(is.na(x)){x=0}else{if(x>0){x=1}else{if(x<0){x=-1}}};x}) %>% as.numeric()) * WMW_sum$WMW_logFDR

write.table(WMW_sum,"../output/ArchR/chromVAR/WMW_res_sgRNA_summary.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#****************************************#
# drawing volcano plot
#****************************************#

for(sgRNA_TMP in sgRNA_list$HTO[3:10]){
  df_tmp        = WMW_sum %>% dplyr::filter(sgRNA==sgRNA_TMP)

  pos_num_tmp = sum((df_tmp$diff_med>0)&(df_tmp$Diff=="significant"))
  neg_num_tmp = sum((df_tmp$diff_med<0)&(df_tmp$Diff=="significant"))

  df_tmp$diff_med[-3>df_tmp$diff_med] = -3
  df_tmp$diff_med[df_tmp$diff_med>3]  = 3
  df_tmp$WMW_logFDR[df_tmp$WMW_logFDR>150]  = 150
  df_tmp = df_tmp %>% arrange(desc(WMW_logFDR))
  df_tmp$TF2 = strsplit(df_tmp$TF,"_") %>% sapply(.,function(x){x[1]})
  TF_list = c("JUN", "FOSL2", "NFATC2", "NR4A2", "NFATC3", "BATF", "NFKB1", "RELA", "ETV2", "ELF1", "IRF2", "ELK1", "ETS1", "RELB", "HIVEP3", "CBFB", "RUNX1", "RUNX2", "RUNX3" ,"NFKB2")
  df_tmp$label = NA
  df_tmp$label[is.element(df_tmp$TF2,TF_list)] = df_tmp$TF[is.element(df_tmp$TF2,TF_list)]

  col_tmp = sgRNA_list$color[sgRNA_list$HTO==sgRNA_TMP]

  p = ggplot(df_tmp,aes(x=diff_med,y=WMW_logFDR,fill=Diff))+
         geom_hline(yintercept=log10(20),linetype="dashed")+
         geom_vline(xintercept=0.25,linetype="dashed")+
         geom_vline(xintercept=-0.25,linetype="dashed")+
         ggrastr::geom_point_rast(size=2,color="black",shape=21)+
         # geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")+
         theme_classic()+
         scale_fill_manual(values=c("grey",col_tmp))+
         scale_x_continuous(breaks=seq(-3,3,by=0.5),limits=c(-3.1,3.1),expand=c(0,0))+
         ## scale_y_continuous(breaks=seq(-50,160,by=50),limits=c(-10,160),expand=c(0,0))+
         theme(axis.text.x  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               axis.text.y  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               axis.title   = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               plot.title   = element_text(colour="black",size=8, face="bold", family = "Helvetica"),
               legend.text  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
               legend.position = "none")+
         labs(title=paste0(sgRNA_TMP," volcano-Plot pos:",pos_num_tmp,", nega:",neg_num_tmp,", full:870 TFs"),x="diff Motif deviation",y="-log10FDR")
  if(is.element(sgRNA_TMP,c("sgCD4_sgGuide1","sgCD4_sgGuide2"))){p = p + scale_y_continuous(breaks=seq(0,60,by=20),limits=c(-5,50),expand=c(0,0))}

  paste0("../output/ArchR/chromVAR/volcanoPlot/",sgRNA_TMP,"_volcanoPlot.pdf") %>% pdf(.,h=5,w=3.8)
   plot(p)
  dev.off()

  p2 = p + geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")

  ############### Fig F / SupFigF ###########################################################################

  paste0("../plots/chromVAR_volcanoPlot/",sgRNA_TMP,"_volcanoPlot_w_repel.pdf") %>% pdf(.,h=5,w=3.8)
   plot(p2)
  dev.off()
}

#****************************************#
# drawing scatter plot
#****************************************#
for(i in 1:4){
  if(i==2){target_sgRNA_1 = sgRNA_list$HTO[2*i +2]
           target_sgRNA_2 = sgRNA_list$HTO[2*i +1]
           }else{
           target_sgRNA_1 = sgRNA_list$HTO[2*i +1]
           target_sgRNA_2 = sgRNA_list$HTO[2*i +2]
          }
  res_tmp   = WMW_sum %>% dplyr::filter(is.element(WMW_sum$sgRNA,c(target_sgRNA_1,target_sgRNA_2)))

  ######
  res_tmp_2 = res_tmp %>% dplyr::select(TF,sgRNA,sign_WMW_logFDR) %>%
               spread( .,key=sgRNA,value=sign_WMW_logFDR)
  res_tmp_2 = res_tmp_2[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_2) = c("TF","target1","target2")
  res_tmp_2$label        = NA
  res_tmp_2              = res_tmp_2 %>% arrange(desc(abs(target1)))
  res_tmp_2$label[1:40] = res_tmp_2$TF[1:40]
  res_tmp_2              = res_tmp_2 %>% arrange(desc(abs(target2)))
  res_tmp_2$label[1:40] = res_tmp_2$TF[1:40]

  ######
  res_tmp_3 = res_tmp %>% dplyr::select(TF,sgRNA,diff_med) %>%
               spread( .,key=sgRNA,value=diff_med)
  res_tmp_3 = res_tmp_3[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_3) = c("TF","target1","target2")
  res_tmp_3$label        = NA
  res_tmp_3              = res_tmp_3 %>% arrange(desc(abs(target1)))
  res_tmp_3$label[1:40] = res_tmp_3$TF[1:40]
  res_tmp_3              = res_tmp_3 %>% arrange(desc(abs(target2)))
  res_tmp_3$label[1:40] = res_tmp_3$TF[1:40]

  ######
  res_tmp_4 = res_tmp %>% dplyr::select(TF,sgRNA,Diff) %>% 
               spread(.,key=sgRNA,value=Diff)
  res_tmp_4 = res_tmp_4[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_4) = c("TF","target1","target2")
  res_tmp_4$diff = "4"
  res_tmp_4$diff[(res_tmp_4$target1=="significant")&    (res_tmp_4$target2=="significant")    ] = "1_target1_target2"
  res_tmp_4$diff[(res_tmp_4$target1=="significant")&    (res_tmp_4$target2=="non_significant")] = "2_target1"
  res_tmp_4$diff[(res_tmp_4$target1=="non_significant")&(res_tmp_4$target2=="significant")    ] = "3_target2"

  ######
  col_tmp = sgRNA_list$color[is.element(sgRNA_list$HTO,c(target_sgRNA_1,target_sgRNA_2))] %>% unique()
  if(length(col_tmp)==1){col          = colorRampPalette(c(col_tmp,"gray"))
                         col_tmp_2    = data.frame(diff=c("1_target1_target2","2_target1","3_target2","4"),col=col(4))
                        }else{
                         col          = colorRampPalette(c(col_tmp[1],col_tmp[2]))
                         col_tmp_2    = data.frame(diff=c("1_target1_target2","2_target1","3_target2","4"),col=c(col(3)[2],col_tmp[1],col_tmp[2],"gray"))
                         }
  col_tmp_list = col_tmp_2[is.element(col_tmp_2$diff,unique(res_tmp_4$diff)),"col"]

  ######
  res_tmp_5 = left_join(res_tmp_2,res_tmp_4 %>% dplyr::select(TF,diff),by="TF")

  p = ggplot(res_tmp_5,aes(x=target1,y=target2,fill=diff))+
         geom_hline(yintercept=0)+
         geom_vline(xintercept=0)+
         ggrastr::geom_point_rast(size=2,color="black",shape=21)+
         theme_minimal()+
         scale_fill_manual(values=c(col_tmp_list))+
         ## scale_x_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         ## scale_y_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         theme(axis.text.x  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.text.y  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.title   = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               plot.title   = element_text(colour="black",size=5, face="bold", family = "Helvetica"),
               legend.text  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               legend.position = "none")+
         labs(title=paste0(target_sgRNA_1," ",target_sgRNA_2," FDR 870 TFs logFDR"),
              x=paste0(target_sgRNA_1, " sign logFDR"),
              y=paste0(target_sgRNA_2, " sign logFDR") )

  paste0("../output/ArchR/chromVAR/scatterPlot_logFDR/",target_sgRNA_1,"_",target_sgRNA_2,"_logFDR_scatterPlot.pdf") %>% pdf(.,h=5,w=5)
   plot(p)
  dev.off()

  ############### SupFig G  ###########################################################################

  p2 = p+geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")+
  paste0("../plots/chromVAR_scatterPlot_logFDR/",target_sgRNA_1,"_",target_sgRNA_2,"_logFDR_scatterPlot_w_label.pdf") %>% pdf(.,h=5,w=5)
   plot(p2)
  dev.off()

  ######
  res_tmp_6 = left_join(res_tmp_3,res_tmp_4 %>% dplyr::select(TF,diff),by="TF")

  p = ggplot(res_tmp_6,aes(x=target1,y=target2,fill=diff))+
         geom_hline(yintercept=0)+
         geom_vline(xintercept=0)+
         ggrastr::geom_point_rast(size=2,color="black",shape=21)+
         theme_minimal()+
         scale_fill_manual(values=c(col_tmp_list))+
         ## scale_x_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         ## scale_y_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         theme(axis.text.x  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.text.y  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.title   = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               plot.title   = element_text(colour="black",size=5, face="bold", family = "Helvetica"),
               legend.text  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               legend.position = "none")+
         labs(title=paste0(target_sgRNA_1," ",target_sgRNA_2," FDR 870 TFs median difference deviation Zscore"),
              x=paste0(target_sgRNA_1, " median difference"),
              y=paste0(target_sgRNA_2, " median difference") )

  paste0("../output/ArchR/chromVAR/scatterPlot_median_diff/",target_sgRNA_1,"_",target_sgRNA_2,"_median_diff_scatterPlot.pdf") %>% pdf(.,h=5,w=5)
   plot(p)
  dev.off()

  ############### SupFig G  ###########################################################################
  p2 = p+geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")+
  paste0("../plots/chromVAR_scatterPlot_median_diff/",target_sgRNA_1,"_",target_sgRNA_2,"_median_diff_scatterPlot_w_label.pdf") %>% pdf(.,h=5,w=5)
   plot(p2)
  dev.off()
}

###################################

target_sgRNA_1 = "sgCD3E_sgGuide2"
target_sgRNA_2 = "sgCD4_sgGuide1"

  res_tmp   = WMW_sum %>% dplyr::filter(is.element(WMW_sum$sgRNA,c(target_sgRNA_1,target_sgRNA_2)))
  
  ######
  res_tmp_2 = res_tmp %>% dplyr::select(TF,sgRNA,sign_WMW_logFDR) %>%
               spread( .,key=sgRNA,value=sign_WMW_logFDR)
  res_tmp_2 = res_tmp_2[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_2) = c("TF","target1","target2")
  res_tmp_2$label        = NA
  res_tmp_2              = res_tmp_2 %>% arrange(desc(abs(target1)))
  res_tmp_2$label[1:40] = res_tmp_2$TF[1:40]
  res_tmp_2              = res_tmp_2 %>% arrange(desc(abs(target2)))
  res_tmp_2$label[1:40] = res_tmp_2$TF[1:40]
  
  ######
  res_tmp_3 = res_tmp %>% dplyr::select(TF,sgRNA,diff_med) %>%
               spread( .,key=sgRNA,value=diff_med)
  res_tmp_3 = res_tmp_3[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_3) = c("TF","target1","target2")
  res_tmp_3$label        = NA
  res_tmp_3              = res_tmp_3 %>% arrange(desc(abs(target1)))
  res_tmp_3$label[1:40] = res_tmp_3$TF[1:40]
  res_tmp_3              = res_tmp_3 %>% arrange(desc(abs(target2)))
  res_tmp_3$label[1:40] = res_tmp_3$TF[1:40]

  ######
  res_tmp_4 = res_tmp %>% dplyr::select(TF,sgRNA,Diff) %>% 
               spread(.,key=sgRNA,value=Diff)
  res_tmp_4 = res_tmp_4[,c("TF",target_sgRNA_1,target_sgRNA_2)]
  colnames(res_tmp_4) = c("TF","target1","target2")
  res_tmp_4$diff = "4"
  res_tmp_4$diff[(res_tmp_4$target1=="significant")&    (res_tmp_4$target2=="significant")    ] = "1_target1_target2"
  res_tmp_4$diff[(res_tmp_4$target1=="significant")&    (res_tmp_4$target2=="non_significant")] = "2_target1"
  res_tmp_4$diff[(res_tmp_4$target1=="non_significant")&(res_tmp_4$target2=="significant")    ] = "3_target2"
  
  ######
  col_tmp_2    = data.frame(diff=c("1_target1_target2","2_target1","3_target2","4"),col=c(sgRNA_list$color[c(6,5,3)],"gray"))
  col_tmp_list = col_tmp_2[is.element(col_tmp_2$diff,unique(res_tmp_4$diff)),"col"]
  
  ######
  res_tmp_5 = left_join(res_tmp_2,res_tmp_4 %>% dplyr::select(TF,diff),by="TF")
  
  p = ggplot(res_tmp_5,aes(x=target1,y=target2,fill=diff))+
         geom_hline(yintercept=0)+
         geom_vline(xintercept=0)+
         ggrastr::geom_point_rast(size=2,color="black",shape=21)+
         theme_minimal()+
         scale_fill_manual(values=c(col_tmp_list))+
         ## scale_x_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         ## scale_y_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         theme(axis.text.x  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.text.y  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.title   = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               plot.title   = element_text(colour="black",size=5, face="bold", family = "Helvetica"),
               legend.text  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               legend.position = "none")+
         labs(title=paste0(target_sgRNA_1," ",target_sgRNA_2," FDR 870 TFs logFDR"),
              x=paste0(target_sgRNA_1, " sign logFDR"),
              y=paste0(target_sgRNA_2, " sign logFDR") )
  
  paste0("../output/ArchR/chromVAR/scatterPlot_logFDR/",target_sgRNA_1,"_",target_sgRNA_2,"_logFDR_scatterPlot.pdf") %>% pdf(.,h=5,w=5)
   plot(p)
  dev.off()
  
  p2 = p+geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")+
  paste0("../output/ArchR/chromVAR/scatterPlot_logFDR/",target_sgRNA_1,"_",target_sgRNA_2,"_logFDR_scatterPlot_w_label.pdf") %>% pdf(.,h=5,w=5)
   plot(p2)
  dev.off()

  ######
  res_tmp_6 = left_join(res_tmp_3,res_tmp_4 %>% dplyr::select(TF,diff),by="TF")

  p = ggplot(res_tmp_6,aes(x=target1,y=target2,fill=diff))+
         geom_hline(yintercept=0)+
         geom_vline(xintercept=0)+
         ggrastr::geom_point_rast(size=2,color="black",shape=21)+
         theme_minimal()+
         scale_fill_manual(values=c(col_tmp_list))+
         ## scale_x_continuous(breaks=seq(-250,200,by=50),limits=c(-260,220),expand=c(0,0))+
         scale_y_continuous(breaks=seq(-2,2,by=1),limits=c(-1.1,1.1),expand=c(0,0))+
         theme(axis.text.x  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.text.y  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               axis.title   = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               plot.title   = element_text(colour="black",size=5, face="bold", family = "Helvetica"),
               legend.text  = element_text(colour="black",size=5,face="bold", family = "Helvetica"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               legend.position = "none")+
         labs(title=paste0(target_sgRNA_1," ",target_sgRNA_2," FDR 870 TFs median difference deviation Zscore"),
              x=paste0(target_sgRNA_1, " median difference"),
              y=paste0(target_sgRNA_2, " median difference") )

  paste0("../output/ArchR/chromVAR/scatterPlot_median_diff/",target_sgRNA_1,"_",target_sgRNA_2,"_median_diff_scatterPlot.pdf") %>% pdf(.,h=5,w=5)
   plot(p)
  dev.off()

  ############### SupFig G  ###########################################################################
  p2 = p+geom_text_repel(aes(label=label), family = "Helvetica",na.rm=TRUE,size=1,col="black")+
  paste0("../plots/chromVAR_scatterPlot_median_diff/",target_sgRNA_1,"_",target_sgRNA_2,"_median_diff_scatterPlot_w_label.pdf") %>% pdf(.,h=5,w=5)
   plot(p2)
  dev.off()

q()
#************************************************************************
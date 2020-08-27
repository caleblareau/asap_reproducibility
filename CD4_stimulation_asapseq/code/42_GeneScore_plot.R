options(stringsAsFactors=F)

#******************************************
#### args        = commandArgs(trailingOnly = T)
#### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP="Perturb_CD4_stim"

#******************************************

suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("ggrepel"))

dir.create("../output/ADT_GeneScore", showWarnings = FALSE)

#******************************************
##### data summarizing
#******************************************

list_bed      = data.frame(PATH=list.files("../output/GeneScore/bed_target","bed$",full.names = TRUE,recursive = TRUE))
list_bed$FILE = basename(list_bed$PATH)
list_bed$cond = list_bed$FILE %>% strsplit(.,"_") %>% sapply(.,function(x){paste(x[1:3],collapse="_")})

for(i in 1:nrow(list_bed)){
    df_tmp = fread(list_bed$PATH[i]) %>% as.data.frame(.)
    colnames(df_tmp) = c("chr","start","end","Gene",list_bed$cond[i])
    if(i==1){target_count=df_tmp
            }else{
             if(!all(target_count$Gene==df_tmp$Gene)){print("Making of target bed file FAILED!!")}
             stopifnot(all.equal(target_count$Gene,df_tmp$Gene))
             target_count = cbind(target_count,df_tmp[,list_bed$cond[i]])
             colnames(target_count)[ncol(target_count)] = list_bed$cond[i]
            }
}

write.table(target_count,"../output/GeneScore/bed_target_fragment_summary.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#####
fragment_size           = read.table("../output/GeneScore/bed_fragment_number_list.txt",header=F,row.names=NULL,stringsAsFactors=F)
colnames(fragment_size) = c("fragment_size","PATH")
fragment_size$colname   = basename(fragment_size$PATH) %>% gsub(".bed$","",.)

#******************************************
##### Normalization -- target gene legion Gene Score --
#******************************************

target_count_f  = target_count[,-c(1:3)] %>% column_to_rownames("Gene") %>% 
                   t(.) %>% as.data.frame(.)
fragment_size_f = data.frame(colname = rownames(target_count_f)) %>%
                   left_join(.,fragment_size,by="colname")

target_count_n =log2( target_count_f / fragment_size_f$fragment_size * 10^4 + 1)

 write.table(target_count_n  %>% rownames_to_column("Name"),"../output/GeneScore/target_GeneScore_Normalized.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#******************************************
##### Normalizing ADT data with CLR
#******************************************
ADT_data_CLR =  fread("../output/Signac/after_filter_Signac/ADT_filtered.txt") %>% 
                 as.data.frame(.) %>% column_to_rownames("CellBarcode")

#******************************************
##### Integration of ADT and GeneScore
#******************************************
corresp_list = fread("../data/ADT_Gene_corresp_list.txt") %>% as.data.frame(.) %>% dplyr::filter(Gene!="")
corresp_list = corresp_list[is.element(corresp_list$Gene,colnames(target_count_n)),]

sgRNA_list   = read.table("../data/sgRNA_list.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")

data_GS = target_count_n %>% rownames_to_column("Name") %>%
           gather(.,key=Gene,value=GS,-Name)
data_GS$HTO         = data_GS$Name %>% strsplit(.,"_") %>% sapply(.,function(x){paste(x[1:2],collapse="_")})
data_GS$CellBarcode = data_GS$Name %>% strsplit(.,"_") %>% sapply(.,function(x){x[3]}) %>% paste0(.,"-1")

data_ADT = ADT_data_CLR[,corresp_list$ADT_name] %>% rownames_to_column("CellBarcode") %>%
            gather(.,key=ADT_name,value=ADT,-CellBarcode)

data_ADT_GS = left_join(data_GS,corresp_list,by="Gene") %>%
               full_join(.,data_ADT,by=c("ADT_name","CellBarcode")) %>%
               dplyr :: select(HTO,CellBarcode,ADT_name,Gene,GS,ADT)

#******************************************
##### Calculation of difference from sgNTC(sgNTCsgGuide1+sgNTCsgGuide2)
#******************************************
data_ADT_GS_sgNTC_stat = data_ADT_GS[grep("sgNTC",data_ADT_GS$HTO),] %>% 
                          group_by(Gene) %>% 
                          summarise_each(funs(mean,sd), GS,ADT) %>% 
                          as.data.frame(.)

ADT_GS_stat_res = c()
for(Gene_tmp in corresp_list$Gene){
  sgNTC_stat_tmp  = data_ADT_GS_sgNTC_stat %>% dplyr::filter(Gene==Gene_tmp)

  data_ADT_GS_tmp = data_ADT_GS %>% dplyr::filter(Gene==Gene_tmp)
  data_ADT_GS_tmp$GS_zscore  = (data_ADT_GS_tmp$GS-sgNTC_stat_tmp$GS_mean)/sgNTC_stat_tmp$GS_sd
  data_ADT_GS_tmp$ADT_zscore = (data_ADT_GS_tmp$ADT-sgNTC_stat_tmp$ADT_mean)/sgNTC_stat_tmp$ADT_sd

  ADT_GS_stat_tmp   = data_ADT_GS_tmp %>%
                       group_by(HTO) %>% 
                       summarise_each(funs(mean),GS_zscore,ADT_zscore) %>% 
                       as.data.frame(.)
  ADT_GS_stat_tmp$Gene = Gene_tmp
  ADT_GS_stat_res = rbind(ADT_GS_stat_res,ADT_GS_stat_tmp)
}

ADT_GS_stat_res = ADT_GS_stat_res[grep("sgNTC",ADT_GS_stat_res$HTO,invert=T),]

write.table(ADT_GS_stat_res,"../output/ADT_GeneScore/ADT_GeneScore_Zscore_w_sgNTC.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

ADT_GS_stat_res$cond = paste0(ADT_GS_stat_res$Gene,"_",ADT_GS_stat_res$HTO) %>% gsub("sgGuide","",.)
ADT_GS_stat_res$label = ""
ADT_GS_stat_res$label[(abs(ADT_GS_stat_res$ADT_zscore)>0.5)|(abs(ADT_GS_stat_res$GS_zscore)>0.2)] = ADT_GS_stat_res$cond[(abs(ADT_GS_stat_res$ADT_zscore)>0.5)|(abs(ADT_GS_stat_res$GS_zscore)>0.2)]

sgRNA_list_wo_sgNTC = sgRNA_list[grep("sgNTC",sgRNA_list$HTO,invert=T),]
ADT_GS_stat_res$HTO = factor(ADT_GS_stat_res$HTO,levels=sgRNA_list_wo_sgNTC$HTO)

p = ggplot(ADT_GS_stat_res,aes(x=ADT_zscore,y=GS_zscore,fill=HTO,label=label))+
     # geom_abline(slope=1,intercept=0,col="gray")+
     geom_vline(xintercept=0,col="black")+
     geom_hline(yintercept=0,col="black")+
     ggrastr::geom_point_rast(size=2,color="black",shape=21)+
     # geom_text_repel(col="black",fontface = "bold",size=1)+
     theme_classic()+
     scale_fill_manual(values=sgRNA_list_wo_sgNTC$color)+
     scale_x_continuous(breaks=seq(-3,2,by=1),limits=c(-2.9,1.9),expand=c(0,0))+
     scale_y_continuous(breaks=seq(-1,1,by=0.2),limits=c(-0.7,0.59),expand=c(0,0))+
     theme(axis.text.x  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
           axis.text.y  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
           axis.title   = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
           plot.title   = element_text(colour="black",size=8, face="bold", family = "Helvetica"),
           legend.text  = element_text(colour="black",size=13,face="bold", family = "Helvetica"),
           legend.position = "none")+
     labs(title=paste0(SAMPLE_TMP," ADT_CLR GeneScore_log w/ sgNTC z-Score Mean"),x="ADT_CLR w/ sgNTC z-Score Mean",y="GeneScore_log w/ sgNTC z-Score Mean")

paste0("../output/ADT_GeneScore/ADT_GeneScore_log_Zscore_w_sgNTC_Mean.pdf") %>% pdf(.,h=5,w=5)
 plot(p)
dev.off()

p2 = p +  geom_text_repel(col="black",fontface = "bold",size=1)

################ Fig G ################################################################################
paste0("../plots/ADT_GeneScore_log_Zscore_w_sgNTC_Mean_w_repel.pdf") %>% pdf(.,h=5,w=5)
 plot(p2)
dev.off()

##########################

ADT_GS_stat_res_2 = ADT_GS_stat_res %>% 
                     dplyr::filter(Gene!="CD3E") %>% dplyr::filter(Gene!="CD4")

print("######### peason correlation between ADT and GeneScore w/o CD3E CD4 #########")
cor.test(ADT_GS_stat_res_2$GS_zscore, ADT_GS_stat_res_2$ADT_zscore, method="pearson") %>% print(.)
####         Pearson's product-moment correlation
#### data:  ADT_GS_stat_res_2$GS_zscore and ADT_GS_stat_res_2$ADT_zscore
#### t = 8.6866, df = 158, p-value = 4.438e-15
#### alternative hypothesis: true correlation is not equal to 0
#### 95 percent confidence interval:
####  0.4533533 0.6650181
#### sample estimates:
####       cor
#### 0.5685213

ADT_GS_stat_res_3 = ADT_GS_stat_res %>% 
                     dplyr::filter(Gene!="CD3E") %>% dplyr::filter(Gene!="CD4") %>% dplyr::filter(Gene!="CD69")

print("######### peason correlation between ADT and GeneScore w/o CD3E CD4 CD69#########")
cor.test(ADT_GS_stat_res_3$GS_zscore, ADT_GS_stat_res_3$ADT_zscore, method="pearson") %>% print(.)
####         Pearson's product-moment correlation
#### data:  ADT_GS_stat_res_3$GS_zscore and ADT_GS_stat_res_3$ADT_zscore
#### t = 11.985, df = 150, p-value < 2.2e-16
#### alternative hypothesis: true correlation is not equal to 0
#### 95 percent confidence interval:
####  0.6079088 0.7725938
#### sample estimates:
####     cor
#### 0.69942


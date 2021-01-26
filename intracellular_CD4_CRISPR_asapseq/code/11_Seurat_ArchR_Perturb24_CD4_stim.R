
options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86)) ## GRCh38
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(BuenColors))
suppressMessages(library(gridExtra))
suppressMessages(library(data.table))
suppressMessages(library(ArchR))

set.seed(111)

#******************************************
args        = commandArgs(trailingOnly = T)
SAMPLE_TMP  = as.character( args[1] )
## SAMPLE_TMP="Perturb24_CD4_stim"

#******************************************
paste0("SAMPLE_TMP : ",SAMPLE_TMP) %>% print()

cellranger_Dir = paste0("../../../asap_large_data_files/intraCD4_stimulation_asapseq/input/")

########################################################
########## Pre-processing workflow            ##########
########################################################

## count
counts   = Read10X_h5(filename = paste0(cellranger_Dir,"filtered_peak_bc_matrix.h5") )
paste0("peak_count ,", " tag : ",ncol(counts),", peak : ",nrow(counts)) %>% print()
## "peak_count , tag : 21,135, peak : 122,662"

## metadata
metadata = read.csv(file = paste0("../data/cellranger/singlecell.csv"),header = TRUE,row.names = 1)
paste0("metadata ,", " tag : ",nrow(metadata)) %>% print()
## "metadata , tag : 1,012,701"

## HTO_ADT data
HTO_ADT_data = read.table("../data/Perturb24_CD4_stim_ADTHTO_summary.txt",header=T,row.names=1,stringsAsFactors=F)
colnames(HTO_ADT_data) = gsub("_","-",colnames(HTO_ADT_data))

HTO_data     = HTO_ADT_data %>% dplyr::select(starts_with("HTO-A"))
surface_data = HTO_ADT_data %>% dplyr::select(starts_with("surface-"))
intra_data   = HTO_ADT_data %>% dplyr::select(starts_with("intra-"))

bc_ATAC      = colnames(counts)
bc_HTO       = rownames(HTO_data)
bc_surface   = rownames(surface_data)
bc_intra     = rownames(intra_data)
metaName     = rownames(metadata)
bc_common    = Reduce(intersect,list(bc_ATAC,bc_HTO,bc_surface,bc_intra,metaName))

#################################
counts_2       = counts[,bc_common]
HTO_data_2     = HTO_data[bc_common,]
surface_data_2 = surface_data[bc_common,]
intra_data_2   = intra_data[bc_common,]
metadata_2     = metadata[bc_common,]

paste0("barcode# : ", length(bc_common))    %>% print() ## "barcode#  :  21,135"
paste0("ATACpeak# : ",nrow(counts_2))       %>% print() ## "ATACpeak# : 122,662"
paste0("Hashtag# : ", ncol(HTO_data_2))     %>% print() ## "Hashtag#  :      14"
paste0("surface# : ", ncol(surface_data_2)) %>% print() ## "surface#  :      53"
paste0("intra# : ",   ncol(intra_data_2))   %>% print() ## "intra#    :      7"

################
fragment.path = paste0(cellranger_Dir,"fragments.tsv.gz")
chrom_assay   = CreateChromatinAssay(counts = counts,
                                     sep = c(":", "-"),
                                     genome = 'hg38',
                                     fragments = fragment.path,
                                     min.cells = 1,
                                     min.features = 200        )
Seurat_obj    = CreateSeuratObject(counts    = chrom_assay,
                                   assay     = "peaks",
                                   meta.data = metadata_2)
print(Seurat_obj)
## An object of class Seurat 
## 122502 features across 21135 samples within 1 assay 
## Active assay: peaks (122502 features, 0 variable features)

#****************

# extract gene annotations from EnsDb
annotations                 = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86) ## EnsDb.Hsapiens.v75
seqlevelsStyle(annotations) = 'UCSC'
genome(annotations)         = "hg38"
Annotation(Seurat_obj)      = annotations

print(Seurat_obj[["peaks"]])
## ChromatinAssay data with 122502 features for 21135 cells
## Variable features: 0 
## Genome: hg38 
## Annotation present: TRUE 
## Motifs present: FALSE 
## Fragment files: 1 

########################################################
########## Computing QC Metrics               ##########
########################################################

## nucleosomal signal
Seurat_obj = NucleosomeSignal(object = Seurat_obj)
Seurat_obj = TSSEnrichment(object = Seurat_obj, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
Seurat_obj$pct_reads_in_peaks = Seurat_obj$peak_region_fragments / Seurat_obj$passed_filters * 100
Seurat_obj$blacklist_ratio    = Seurat_obj$blacklist_region_fragments / Seurat_obj$peak_region_fragments

## TSS enrichment
Seurat_obj$high.tss = ifelse(Seurat_obj$TSS.enrichment > 2, 'High', 'Low')

paste0("../output/Signac/QC_plot/",SAMPLE_TMP,"_TSSplot.pdf") %>%  pdf(.,h=5,w=10)
 TSSPlot(Seurat_obj, group.by = 'high.tss') + ggtitle(paste0(SAMPLE_TMP," TSS enrichment score")) + NoLegend()
dev.off()

## Fragment enrichment
Seurat_obj$nucleosome_group <- ifelse(Seurat_obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

paste0("../output/Signac/QC_plot/",SAMPLE_TMP,"_nucleosome_group_plot.pdf") %>% pdf(.,h=5,w=10)
 FragmentHistogram(object = Seurat_obj, group.by = 'nucleosome_group')+ plot_annotation( title = SAMPLE_TMP ) + NoLegend()
dev.off()

## QC summary
p1 = VlnPlot(Seurat_obj, c('pct_reads_in_peaks', 'peak_region_fragments',
                           'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'), pt.size = 0, ncol = 5)
## p2 <- VlnPlot(Seurat_obj, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()

pdf(paste0("../output/Signac/QC_plot/",SAMPLE_TMP,"_QC_plot.pdf"),h=5,w=12.5)
 p1 + plot_annotation( title = SAMPLE_TMP)
dev.off()

###################
### ###################
# Normalize HTO/ADT data, here we use centered log-ratio (CLR) transformation
Seurat_obj[["HTO"]] <- CreateAssayObject(counts = t(HTO_data_2))
Seurat_obj          <- NormalizeData(Seurat_obj, assay = "HTO", normalization.method = "CLR")
Seurat_obj          <- ScaleData( Seurat_obj, assay = "HTO")

### ###################
Seurat_obj[["ADTsurface"]] <- CreateAssayObject(counts = t(surface_data_2))
Seurat_obj                 <- NormalizeData(Seurat_obj, assay = "ADTsurface", normalization.method = "CLR")
Seurat_obj                 <- ScaleData(Seurat_obj, assay = "ADTsurface")

### ###################
Seurat_obj[["ADTintra"]] <- CreateAssayObject(counts = t(intra_data_2))
Seurat_obj               <- NormalizeData(Seurat_obj, assay = "ADTintra", normalization.method = "CLR")
Seurat_obj               <- ScaleData(Seurat_obj, assay = "ADTintra")

###################
######### Filter out noise information.
Seurat_obj_pass <- subset(x = Seurat_obj,
                          subset = peak_region_fragments < 35000 &
                                   passed_filters < 60000 & ## 20000
                                   passed_filters > 1000 &
                                   pct_reads_in_peaks > 30 &
                                   blacklist_ratio < 0.05 &
                                   nucleosome_signal < 2.5 & ## 4
                                   TSS.enrichment > 2 )

print(Seurat_obj_pass)
## An object of class Seurat 
## 122576 features across 19766 samples within 4 assays 
## Active assay: peaks (122502 features, 0 variable features)
##  3 other assays present: HTO, ADTsurface, ADTintra

cell_barcode_list = colnames(Seurat_obj_pass)

# omit outliers which have too much ADT read counts. ############################# therhold == 0.999
positive_ADTsurface_tmp = GetAssayData(Seurat_obj_pass, assay = "ADTsurface", slot = 'counts') %>% as.matrix()
positive_ADTintra_tmp   = GetAssayData(Seurat_obj_pass, assay = "ADTintra",   slot = 'counts') %>% as.matrix()
positive_ADT_tmp = rbind(positive_ADTsurface_tmp,positive_ADTintra_tmp)

taget_ctrl_ADT = c("surface-mouseIgG1","intra-mouseIgG2")

outlier_list = c()
for(taget_ctrl_tmp in taget_ctrl_ADT){
  threshold_tmp = quantile(positive_ADT_tmp[taget_ctrl_tmp,],.99) %>% round()         ## .999 .998 .997 .996 .995
  outlier_tmp = colnames(positive_ADTsurface_tmp)[positive_ADT_tmp[taget_ctrl_tmp,]>=threshold_tmp]
  outlier_list = c(outlier_list,outlier_tmp)
}

outlier_list_sum = outlier_list %>% table(.)
outlier_ADT = names(outlier_list_sum)[outlier_list_sum>=2]
### 28 cell barcodes

data.frame(ADT_outlier=outlier_ADT) %>%
 write.table(.,paste0("../output/Signac/after_filter_Signac/ADT_isotype_outlier.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)
Seurat_obj_pass_2  = Seurat_obj_pass[,!is.element(colnames(Seurat_obj_pass),outlier_ADT)]

################################################

Seurat_obj_pass_2_p1 = subset(Seurat_obj_pass_2, cells=grep("-1",colnames(Seurat_obj_pass_2), value=T))
Seurat_obj_pass_2_p2 = subset(Seurat_obj_pass_2, cells=grep("-2",colnames(Seurat_obj_pass_2), value=T))

#######################
Seurat_obj_pass_2_p1 = NormalizeData(Seurat_obj_pass_2_p1, assay = "HTO", normalization.method = "CLR")
Seurat_obj_pass_2_p1 = ScaleData(Seurat_obj_pass_2_p1, assay = "HTO")
Seurat_obj_pass_2_p1 = HTODemux(Seurat_obj_pass_2_p1,  assay = "HTO", positive.quantile = 0.99) # .999

table(Seurat_obj_pass_2_p1$HTO_classification.global)
## Doublet Negative  Singlet
##    8866      21     601 

#######################
Seurat_obj_pass_2_p2 = NormalizeData(Seurat_obj_pass_2_p2, assay = "HTO", normalization.method = "CLR")
Seurat_obj_pass_2_p2 = ScaleData(Seurat_obj_pass_2_p2, assay = "HTO")
Seurat_obj_pass_2_p2 = HTODemux(Seurat_obj_pass_2_p2,  assay = "HTO", positive.quantile = 0.99) # .999

table(Seurat_obj_pass_2_p2$HTO_classification.global)
## Doublet Negative  Singlet 
##    9616       49      585

#######################

HTO_class_data = rbind(Seurat_obj_pass_2_p1@meta.data %>% dplyr::select("HTO_maxID","HTO_secondID","HTO_margin","HTO_classification","HTO_classification.global","hash.ID"),
                       Seurat_obj_pass_2_p2@meta.data %>% dplyr::select("HTO_maxID","HTO_secondID","HTO_margin","HTO_classification","HTO_classification.global","hash.ID"))
HTO_class_data = HTO_class_data[colnames(Seurat_obj_pass_2),]

Seurat_obj_pass_2$HTO_maxID                 = HTO_class_data$HTO_maxID
Seurat_obj_pass_2$HTO_secondID              = HTO_class_data$HTO_secondID
Seurat_obj_pass_2$HTO_margin                = HTO_class_data$HTO_margin
Seurat_obj_pass_2$HTO_classification        = HTO_class_data$HTO_classification
Seurat_obj_pass_2$HTO_classification.global = HTO_class_data$HTO_classification.global
Seurat_obj_pass_2$hash.ID                   = HTO_class_data$hash.ID

###################
HTO_convert = data.frame(CellBarcode=names(Seurat_obj_pass_2$HTO_classification),
                         HTO_classification=Seurat_obj_pass_2$HTO_classification)
sgRNA_list = read.csv("../data/sgRNA_HTO_list.csv",header=T,row.names=NULL,stringsAsFactors=F)
sgRNA_list_2 = sgRNA_list %>% dplyr::select("HTO_2","color") %>% unique()

HTO_convert = left_join(HTO_convert,sgRNA_list,by="HTO_classification")
HTO_convert$HTO[is.na(HTO_convert$HTO)] = "Negative"
HTO_convert$HTO_2[is.na(HTO_convert$HTO_2)] = "Negative"

write.table(HTO_convert,"../output/Signac/after_filter_Signac/HTO_convert_list.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

### reading color list
Seurat_obj_pass_2$HTO_final   = factor(HTO_convert$HTO,levels=c(sgRNA_list$HTO,"Negative"))
Seurat_obj_pass_2$HTO_final_2 = factor(HTO_convert$HTO_2,levels=c(sgRNA_list_2$HTO_2,"Negative"))

ASAP = subset(x = Seurat_obj_pass_2,
                      subset = HTO_final_2!="Negative" )
ASAP$HTO_final     = factor(ASAP$HTO_final,levels=sgRNA_list$HTO       )
ASAP$HTO_final_v2  = factor(ASAP$HTO_final_2,levels=sgRNA_list_2$HTO_2 )

########################################################
########## Qunatile normalization of ADT
########################################################

ADTintra_data   = GetAssayData(ASAP, assay = "ADTintra",   slot = "data") %>% t() %>% as.data.frame()
ADTsurface_data = GetAssayData(ASAP, assay = "ADTsurface", slot = "data") %>% t() %>% as.data.frame()
ADT_data = cbind(ADTintra_data,ADTsurface_data)

tmp_list = data.frame(CellBarcode_2=names(ASAP$HTO_final_2),HTO_final_2=ASAP$HTO_final_2)

for(iii in 1:ncol(ADT_data)){
  for(tmp_HTO in sgRNA_list_2$HTO){
    tmp_CellBarcode = tmp_list$CellBarcode_2[tmp_list$HTO_final_2==tmp_HTO]
    tmp_result = ADT_data[is.element(rownames(ADT_data),tmp_CellBarcode),iii]
    tmp_thresh = quantile(tmp_result,.99)
    tmp_result[tmp_result>tmp_thresh] = tmp_thresh
    ADT_data[is.element(rownames(ADT_data),tmp_CellBarcode),iii] = tmp_result
  }
}

ASAP[["ADTintraQuan"]]   <- CreateAssayObject(counts = ADT_data %>% dplyr::select(starts_with("intra")) %>% t())
ASAP[["ADTsurfaceQuan"]] <- CreateAssayObject(counts = ADT_data %>% dplyr::select(starts_with("surface")) %>% t())

########################################################
########## VlnPlot of HTO / ADT
########################################################

################
Idents(ASAP) <- "HTO_final"

##
for(tmp_num in 1:nrow(ASAP[["ADTsurfaceQuan"]])){
    target_tmp = rownames(ASAP[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=sgRNA_list$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_gRNA/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in 1:nrow(ASAP[["ADTintraQuan"]])){
    target_tmp = rownames(ASAP[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTintraQuan", features = target_tmp, ,slot = "data",same.y.lims = TRUE,cols=sgRNA_list$color, pt.size = 0)+ NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_gRNA/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

################################################################################################################## plot of Figure7f
sgRNA_list_tmp = sgRNA_list[c(1:4,7,8,15,16,13,14),]
Seurat_obj_tmp = ASAP[,is.element(ASAP$HTO_final,sgRNA_list_tmp$HTO)]
Seurat_obj_tmp$HTO_final  = factor(Seurat_obj_tmp$HTO_final,levels=sgRNA_list_tmp$HTO       )
Idents(Seurat_obj_tmp) <- "HTO_final"

for(tmp_num in c(7,28)){
    target_tmp = rownames(Seurat_obj_tmp[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(Seurat_obj_tmp, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=sgRNA_list_tmp$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/Figure7f/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=4,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in c(1,6)){
    target_tmp = rownames(Seurat_obj_tmp[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(Seurat_obj_tmp, assay = "ADTintraQuan", features = target_tmp, ,slot = "data",same.y.lims = TRUE,cols=sgRNA_list_tmp$color, pt.size = 0)+ NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/Figure7f/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=4,w=7)
      plot(p1)
    dev.off()
}


########################################################
########## VlnPlot of HTO / ADT
########################################################

Idents(ASAP) <- "HTO_final_v2"

##
for(tmp_num in 1:nrow(ASAP[["ADTsurfaceQuan"]])){
    target_tmp = rownames(ASAP[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=sgRNA_list_2$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_gene/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in 1:nrow(ASAP[["ADTintraQuan"]])){
    target_tmp = rownames(ASAP[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTintraQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=sgRNA_list_2$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_gene/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

################################################################################################################## plot of FigureSupple7d
for(tmp_num in c(26,36,43)){
    target_tmp = rownames(ASAP[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=sgRNA_list_2$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/FigureSupple7d/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in 4){
    target_tmp = rownames(ASAP[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTintraQuan", features = target_tmp, ,slot = "data",same.y.lims = TRUE,cols=sgRNA_list_2$color, pt.size = 0)+ NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/FigureSupple7d/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

#**********************************************
## HTO heatmap
#**********************************************
################################################################################################################## plot of FigureSupple7c

HTO_data = GetAssayData(ASAP, assay = "HTO",  slot = "data")       %>% t() %>% as.data.frame()

tmp_thresh = c(strsplit(sgRNA_list$HTO_classification,"_") %>% sapply(.,function(x){x[1]}),strsplit(sgRNA_list$HTO_classification,"_") %>% sapply(.,function(x){x[2]})) %>% table() %>% as.data.frame()
colnames(tmp_thresh)[1]="HTO"

tmp_thresh$thresh = ((1 - tmp_thresh$Freq/nrow(sgRNA_list))*100) %>% ceiling(.)/100
tmp_thresh$thresh = (1+tmp_thresh$thresh)/2
tmp_thresh_2 = left_join(data.frame(HTO=colnames(HTO_data)),tmp_thresh,by="HTO")

for(tmp_num in 1:ncol(HTO_data)){
    HTO_data[,tmp_num] = 4*HTO_data[,tmp_num]/quantile(HTO_data[,tmp_num],tmp_thresh_2$thresh[tmp_num])
}

HTO_data[HTO_data<1.5] = 0
HTO_data[HTO_data>4] = 4

hclust_order = function(DATA,METHOD="ward.D2"){
  suppressMessages(library(ggdendro))
  rd     = dist(DATA)
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  col_order = as.character(ddata$labels$label)
  col_order
}

####################################################################################### 
################################################################################################################## plot of FigureSupple7c

tmp_list = list()
for(tmp_num in 1:nrow(sgRNA_list)){
    tmp_list[[tmp_num]] = hclust_order(HTO_data[colnames(ASAP)[ASAP$HTO_final==sgRNA_list$HTO[tmp_num]],])
}
name_list = unlist(tmp_list)

data = HTO_data[name_list,] %>% t() %>% as.data.frame()

col_list_HTOheatmap = colorRampPalette(brewer.pal(9,"Blues"))
df = data %>% rownames_to_column("rownames") %>% gather(.,key=colnames,value=score,-rownames)

df2 = left_join(data.frame(CellBarcode=colnames(data)),ASAP@meta.data %>% dplyr::select("HTO_final") %>% rownames_to_column("CellBarcode"),by="CellBarcode")
df2$HTO_final = factor(df2$HTO_final,levels=sgRNA_list$HTO)

p1  = ggplot(df,aes(x=colnames,y=rownames,fill=score))+
    ggrastr::geom_tile_rast()+
     scale_fill_gradientn(colours=col_list_HTOheatmap(250))+ 
    theme_bw() +
    scale_x_discrete(limits=unique(df$colnames)) +
    scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
          axis.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
          strip.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="Cell Barcode",y="HTO")

p2 = ggplot(df2, aes(x=CellBarcode,y=1,fill=HTO_final))+
      ggrastr::geom_tile_rast()+
      scale_fill_manual(values = sgRNA_list$color)+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_discrete(limits=colnames(data)) +
      theme(axis.title.x=element_blank(),axis.ticks=element_blank(),
            axis.text=element_blank(),
            axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
            plot.title = element_text(size=10,face="bold", family = "Helvetica"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
      labs(title=paste0(SAMPLE_TMP," HTO"),y="sgRNA")

paste0("../plots/FigureSupple7c/",SAMPLE_TMP,"_HTO_heatmap.pdf") %>% pdf(.,h=5,w=15)
 p2/p1+plot_layout(ncol=1, heights=c(1,9))
dev.off()

########################################################################################################################

####################################################################
##############  UMAP with peak data 
####################################################################
################################################################################################################## plot of Figure7e

DefaultAssay(ASAP) <- "peaks"

#Dimensional Reduction
ASAP <- RunTFIDF(ASAP)
ASAP <- FindTopFeatures(ASAP, min.cutoff = 'q0')
ASAP <- RunSVD(ASAP)

#UMAP
ASAP  <- RunUMAP(object = ASAP, reduction = 'lsi',  reduction.name = 'peak.umap', reduction.key = 'peakUMAP_',dims = 2:30,seed.use=1)
ASAP  <- FindNeighbors(object = ASAP, reduction = 'lsi', dims = 2:30)
ASAP  <- FindClusters(object = ASAP, verbose = FALSE, algorithm = 3)

cluster_list = read.csv("../data/UMAP_cluster_list.csv",header=T,row.names=NULL,stringsAsFactors=F)
ASAP$FinalCluster = ASAP$seurat_clusters %>%
                       gsub("^0$","FinalC3",.) %>%
                       gsub("^1$","FinalC5",.) %>%
                       gsub("^2$","FinalC2",.) %>%
                       gsub("^3$","FinalC4",.) %>%
                       gsub("^4$","FinalC1",.) %>%
                       gsub("^5$","FinalC6",.) %>%
                       gsub("^6$","FinalC3",.) %>%
                       factor(.,levels=cluster_list$FinalCluster)

p1    <- DimPlot(object = ASAP, reduction="peak.umap", label = TRUE, group.by="FinalCluster",cols=cluster_list$color) + NoLegend() + ggtitle("ASAP_peak_UMAP cluster")

pdf(paste0("../plots/Figure7e/",SAMPLE_TMP,"_UMAP_ASAPpeaks.pdf"), w = 5, h = 5)
 plot(p1)
dev.off()

#############################################

meta_data  = ASAP@meta.data
write.table(meta_data %>% rownames_to_column("CellBarcode") %>% dplyr::select(c("CellBarcode","FinalCluster","HTO_final","HTO_final_2")),
            "../output/Signac/after_filter_Signac/Seurat_cluster_data.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

table_cluster = table(meta_data$HTO_final,meta_data$FinalCluster) %>% prop.table(.,1) %>%
                as.data.frame()
colnames(table_cluster) = c("HTO","FinalCluster","Freq")

table_cluster$HTO          = factor(table_cluster$HTO,          levels = sgRNA_list$HTO            %>% rev() )
table_cluster$FinalCluster = factor(table_cluster$FinalCluster, levels = cluster_list$FinalCluster %>% rev() )

p1 = ggplot(table_cluster ,aes(x=HTO,y=Freq,fill=FinalCluster))+
        geom_bar(stat = "identity")+theme_bw()+
        coord_flip()+
        scale_fill_manual(values = cluster_list$color %>% rev())+
        theme(axis.text   = element_text(size=10),
              axis.title  = element_blank(),
              axis.text.y = element_text(size=10,face="bold",colour="black", family = "Helvetica"),
              axis.text.x = element_text(size=10,hjust=1,vjust=1,face="bold",colour="black", family = "Helvetica"),
              legend.position="none")+
        labs(title=paste0(SAMPLE_TMP," Cluster proportion"))

paste0("../plots/Figure7e/",SAMPLE_TMP,"_each_cluster_proportion.pdf") %>% pdf(.,h=5,w=5)
 plot(p1)
dev.off()

#####

dimplot_color_list   = jdb_palette("solar_extra",type="continuous") %>% as.character(.)

DefaultAssay(ASAP) = "ADTsurface"
for(tmp_num in 1:nrow(ASAP[["ADTsurface"]])){
  target_tmp = rownames(ASAP[["ADTsurface"]])[tmp_num]
  p1 = FeaturePlot(ASAP, features = target_tmp , min.cutoff = "q10",  max.cutoff = "q90",
                  cols = dimplot_color_list, pt.size=0.1, order = TRUE, reduction = "peak.umap")
  p1 = AugmentPlot(plot=p1,width=4,height=4,dpi=100)

  paste0("../output/Signac/DimPlot/",SAMPLE_TMP,"_ADTsurface_DimPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>% pdf(., w = 5, h = 5)
   plot(p1)
  dev.off()
}

DefaultAssay(ASAP) = "ADTintra"
for(tmp_num in 1:nrow(ASAP[["ADTintra"]])){
  target_tmp = rownames(ASAP[["ADTintra"]])[tmp_num]
  p1 = FeaturePlot(ASAP, features = target_tmp , min.cutoff = "q10",  max.cutoff = "q90",
                  cols = dimplot_color_list, pt.size=0.1, order = TRUE, reduction = "peak.umap")
  p1 = AugmentPlot(plot=p1,width=4,height=4,dpi=100)

  paste0("output/Signac/DimPlot/",SAMPLE_TMP,"_ADTintra_DimPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>% pdf(., w = 5, h = 5)
   plot(p1)
  dev.off()
}

####################
################################################################################################################## plot of Figure7g

DefaultAssay(ASAP) = "ADTsurface"
for(tmp_num in c(7,28) ){
  target_tmp = rownames(ASAP[["ADTsurface"]])[tmp_num]
  p1 = FeaturePlot(ASAP, features = target_tmp , min.cutoff = "q10",  max.cutoff = "q90",
                  cols = dimplot_color_list, pt.size=0.1, order = TRUE, reduction = "peak.umap")
  p1 = AugmentPlot(plot=p1,width=4,height=4,dpi=100)

  paste0("../plots/Figure7g/",SAMPLE_TMP,"_ADTsurface_DimPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>% pdf(., w = 5, h = 5)
   plot(p1)
  dev.off()
}

DefaultAssay(ASAP) = "ADTintra"
for(tmp_num in c(1,6) ){
  target_tmp = rownames(ASAP[["ADTintra"]])[tmp_num]
  p1 = FeaturePlot(ASAP, features = target_tmp , min.cutoff = "q10",  max.cutoff = "q90",
                  cols = dimplot_color_list, pt.size=0.1, order = TRUE, reduction = "peak.umap")
  p1 = AugmentPlot(plot=p1,width=4,height=4,dpi=100)

  paste0("../plots/Figure7g/",SAMPLE_TMP,"_ADTintra_DimPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>% pdf(., w = 5, h = 5)
   plot(p1)
  dev.off()
}

########################################################
########## VlnPlot of HTO / ADT
########################################################

Idents(ASAP) <- "FinalCluster"

##
for(tmp_num in 1:nrow(ASAP[["ADTsurfaceQuan"]])){
    target_tmp = rownames(ASAP[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=cluster_list$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_cluster/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in 1:nrow(ASAP[["ADTintraQuan"]])){
    target_tmp = rownames(ASAP[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTintraQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=cluster_list$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../output/Signac/VlnPlot/each_cluster/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

################################################################################################################## plot of Figure7h

for(tmp_num in c(4,19,31,43) ){
    target_tmp = rownames(ASAP[["ADTsurfaceQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTsurfaceQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=cluster_list$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/Figure7h/",SAMPLE_TMP,"_ADTsurface_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

for(tmp_num in c(4,6) ){
    target_tmp = rownames(ASAP[["ADTintraQuan"]])[tmp_num]
    p1 = VlnPlot(ASAP, assay = "ADTintraQuan", features = target_tmp, slot = "data",same.y.lims = TRUE,cols=cluster_list$color, pt.size = 0) + NoLegend() + plot_annotation( title = paste0(SAMPLE_TMP," ",target_tmp))
    paste0("../plots/Figure7h/",SAMPLE_TMP,"_ADTintra_VlnPlot_",formatC(tmp_num,width=2,flag="0"),"_",target_tmp, ".pdf") %>%  pdf(.,h=5,w=7)
      plot(p1)
    dev.off()
}

saveRDS(object=ASAP, file= paste0("../output/Signac/after_filter_Signac/",SAMPLE_TMP,"_Seurat_object_tmp.rds"))

########################################################
########## ArchR
########################################################

rownames(meta_data) = paste0(SAMPLE_TMP,"#",rownames(meta_data))

addArchRThreads(threads = 10)
addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles      = fragment.path ,
  sampleNames     = SAMPLE_TMP ,
  filterTSS       = 1, #Dont set this too high because you can always increase later
  filterFrags     = 500, 
  addTileMat      = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE,
  QCDir   = paste0("../output/ArchR/QC"), 
  logFile = paste0("../output/ArchR/log/",SAMPLE_TMP,"_createArrowFiles.log")
)

doubScores <- addDoubletScores(input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  outDir  = paste0("../output/ArchR/QC"),
  logFile = paste0("../output/ArchR/log/",SAMPLE_TMP,"_addDoubletScores.log")
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0("../output/ArchR/",SAMPLE_TMP),
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

## cell barcode after QC filter by Seurat
target_cell = intersect(rownames(meta_data),rownames(proj))
proj_subset = subsetCells(ArchRProj = proj, cellNames = target_cell)

GeneScoreMat = getMatrixFromProject(ArchRProj = proj_subset,
                                    useMatrix = "GeneScoreMatrix",
                                    useSeqnames = NULL,
                                    verbose = TRUE,
                                    binarize = FALSE,
                                    threads = getArchRThreads(),
                                    logFile = paste0("output/ArchR/log/",SAMPLE_TMP,"_getMatrixFromProject.log") )

data_GeneScore = assays(GeneScoreMat)$GeneScoreMatrix %>% as.data.frame()
rownames(data_GeneScore) = rowData(GeneScoreMat)$name
colnames(data_GeneScore) = gsub(paste0(SAMPLE_TMP,"#"),"",colnames(data_GeneScore))
data_GeneScore_2 = data_GeneScore[,colnames(ASAP)]

ASAP[["GeneScore"]] = CreateAssayObject(counts = data_GeneScore_2 )
ASAP                = NormalizeData(ASAP, assay="GeneScore",normalization.method = "LogNormalize", scale.factor = 10000)

#############################  omit outliers

data_GeneScore   = GetAssayData(ASAP, assay = "GeneScore",   slot = "data") %>% as.data.frame() %>% t() 
## data_GeneScore   = data_GeneScore[,c("IFNG","GZMB","MKI67")] ##################################################33

for(iii in 1:ncol(data_GeneScore)){
    tmp_result = data_GeneScore[,iii]
    tmp_thresh_99 = quantile(tmp_result,.99)
    tmp_result[tmp_result>tmp_thresh_99] = tmp_thresh_99
    data_GeneScore[,iii] = tmp_result
}

################################################################################################################## plot of FigureSupple7e
ASAP[["GeneScoreQuan"]]   <- CreateAssayObject(data = data_GeneScore %>% t())

Idents(ASAP) <- "FinalCluster"
for(target_tmp in c("IFNG","GZMB")){
  p1 = VlnPlot(ASAP, assay = "GeneScoreQuan", features = target_tmp,slot = "data",same.y.lims = FALSE, pt.size = 0,cols=cluster_list$color) + NoLegend() +plot_annotation( title = paste0(SAMPLE_TMP," GeneScore"))
  paste0("../plots/FigureSupple7e/",SAMPLE_TMP,"_GeneScore_VlnPlot_",target_tmp, ".pdf") %>%  pdf(.,h=3,w=5)
    plot(p1)
  dev.off()
}

########################################################
### Saving Data After Filter                         ###
########################################################
saveRDS(object=ASAP, file= paste0("output/Signac/after_filter_Signac/",SAMPLE_TMP,"_Seurat_object_final.rds"))

########################################################
### Saving sessionInfo                               ###
########################################################
sink(paste0("../output/sessionInfo.txt"))
sessionInfo()
sink()






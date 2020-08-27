options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(BUSpaRse))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86)) ## GRCh38
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(dsb))
suppressMessages(library(magrittr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(BuenColors))
suppressMessages(library(gridExtra))

set.seed(111)

#******************************************
SAMPLE_TMP     = "Perturb_CD4_stim"
HashTag_list   = c("sgNTC","sgCD4","sgCD3E","sgZAP70","sgNFKB2","sgGuide1","sgGuide2")

#dir_list = c("output/Signac/after_DSB","output/Signac/after_filter_Signac","output/Signac/QC_plot",
#             "output/Signac/RidgePlot","output/Signac/VlnPlot","plots")
#for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}

########################################################
########## Pre-processing workflow            ##########
########################################################

## count
counts   = Read10X_h5(filename = paste0("../data/filtered_peak_bc_matrix.h5") )
paste0("peak_count ,", " tag : ",ncol(counts),", peak : ",nrow(counts)) %>% print()
#### "peak_count , tag : 9151, peak : 94838"

## metadata
metadata = read.csv(file = paste0("../data/singlecell.csv"),header = TRUE,row.names = 1)
paste0("metadata ,", " tag : ",nrow(metadata)) %>% print()
#### "metadata , tag : 549584"

#*************************************************

ADT_data = read_count_output("../data/kallisto/adt", name = "featurecounts", tcc = FALSE) %>% 
           t() %>% as.data.frame()
HTO_data = read_count_output("../data/kallisto/hto", name = "featurecounts", tcc = FALSE) %>%  
          t() %>% as.data.frame()

hashtag_list = read.table("../data/kallisto/hashtag_list.txt",header=T,row.names=NULL,stringsAsFactors=F)

rownames(ADT_data) = paste0(rownames(ADT_data),"-1")
colnames(ADT_data) = gsub("_","--",colnames(ADT_data)) %>% strsplit(.,"--") %>% 
                      sapply(.,function(x){paste(x[2:3],collapse="--")}) %>% gsub("--NA","",.) %>% gsub("--","_",.) 

HTO_data  = HTO_data[,hashtag_list$hashtag]
colnames(HTO_data) = hashtag_list$sgRNA
rownames(HTO_data) = paste0(rownames(HTO_data),"-1")

paste0("HashTag_check : ", all.equal(colnames(HTO_data)%>%sort(),HashTag_list%>%sort()))
#### "HashTag_check : TRUE"

########################################################
########## cell barcoding check               ##########
########################################################

barcode_ATAC  = colnames(counts)
barcode_HTO   = rownames(HTO_data)
barcode_ADT   = rownames(ADT_data)
metaName = rownames(metadata)
barcode_common = Reduce(intersect,list(barcode_ATAC,barcode_HTO,barcode_ADT,metaName))

#################################
counts_2    = counts[,barcode_common]
HTO_data_2  = HTO_data[barcode_common,]
ADT_data_2  = ADT_data[barcode_common,]
metadata_2  = metadata[barcode_common,]

paste0("barcode# : ", length(barcode_common)) %>% print() ## "barcode# : 9151"
paste0("ATACpeak# : ",nrow(counts_2))         %>% print() ## "ATACpeak# : 94838"
paste0("Hashtag# : ", ncol(HTO_data_2))       %>% print() ## "Hashtag# : 7"
paste0("ADT# : ",     ncol(ADT_data_2))       %>% print() ## "ADT# : 37"

########################################################
########## making of Seurat Object            ##########
########################################################

Seurat_obj    = CreateSeuratObject(counts = counts_2,assay = 'peaks',
                                    project = 'ATAC',min.cells = 1, meta.data = metadata_2)
fragment.path = paste0("../../../asap_large_data_files/CD4_stimulation_asapseq/cellranger/fragments.tsv.gz")
Seurat_obj    = SetFragments(object = Seurat_obj,file = fragment.path)
print(Seurat_obj)
#### An object of class Seurat
#### 94523 features across 9151 samples within 1 assay
#### Active assay: peaks (94523 features, 0 variable features)


########################################################
########## Computing QC Metrics               ##########
########################################################

## nucleosomal signal
Seurat_obj <- NucleosomeSignal(object = Seurat_obj)
Seurat_obj$pct_reads_in_peaks <- Seurat_obj$peak_region_fragments / Seurat_obj$passed_filters * 100
Seurat_obj$blacklist_ratio    <- Seurat_obj$blacklist_region_fragments / Seurat_obj$peak_region_fragments

### If we used scale_y_log10(), blacklist_region_fragments=0 became NaNs...
### Seurat_obj_tmp <- subset(x = Seurat_obj,
###                          subset = blacklist_region_fragments >0 )

p1 <- VlnPlot(Seurat_obj, c('pct_reads_in_peaks', 'peak_region_fragments'), pt.size = 0.1)
p2 <- VlnPlot(Seurat_obj, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()

pdf("../output/Signac/QC_plot/QC_plot.pdf",h=5,w=10)
 (p1 | p2) + plot_annotation( title = SAMPLE_TMP)
dev.off()

###################
# create granges object with TSS positions
gene.ranges <- genes(EnsDb.Hsapiens.v86) ###########GRCh38
seqlevelsStyle(gene.ranges) <- 'UCSC'
gene.ranges <- gene.ranges[gene.ranges$gene_biotype == 'protein_coding', ]
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

seqlevelsStyle(tss.ranges) <- 'UCSC'
tss.ranges <- keepStandardChromosomes(tss.ranges, pruning.mode = 'coarse')

# TSS : top 2000 Genes
Seurat_obj <- TSSEnrichment(object = Seurat_obj, tss.positions = tss.ranges[1:2000,])

###################
Seurat_obj$high.tss <- ifelse(Seurat_obj$TSS.enrichment > 2, 'High', 'Low')

pdf("../output/Signac/QC_plot/TSSplot.pdf",h=5,w=10)
 TSSPlot(Seurat_obj, group.by = 'high.tss') + ggtitle(paste0(SAMPLE_TMP," TSS enrichment score")) + NoLegend()
dev.off()

########################################################
########## detecting cells                    ##########
########################################################

# Normalize HTO/ADT data, here we use centered log-ratio (CLR) transformation
Seurat_obj[["HTO"]] <- CreateAssayObject(counts = t(HTO_data_2)[HashTag_list,])
Seurat_obj          <- NormalizeData(Seurat_obj, assay = "HTO", normalization.method = "CLR")
Seurat_obj          <- ScaleData(Seurat_obj, assay = "HTO")
Seurat_obj          <- HTODemux(Seurat_obj, assay = "HTO", positive.quantile = 0.9)

### ###################
Seurat_obj[["ADT"]] <- CreateAssayObject(counts = t(ADT_data_2))
Seurat_obj          <- NormalizeData(Seurat_obj, assay = "ADT", normalization.method = "CLR")
Seurat_obj          <- ScaleData(Seurat_obj, assay = "ADT")

###################
######### Filter out noise information.
Seurat_obj_pass <- subset(x = Seurat_obj,
                          subset = peak_region_fragments > 2000 &
                                   peak_region_fragments < 35000 &
                                   passed_filters < 60000 &
                                   pct_reads_in_peaks > 15 &
                                   blacklist_ratio < 0.05 &
                                   nucleosome_signal < 10 &
                                   TSS.enrichment > 2 )
print(Seurat_obj_pass)
#### An object of class Seurat
#### 94567 features across 6806 samples within 3 assays
#### Active assay: peaks (94523 features, 0 variable features)
####  2 other assays present: HTO, ADT

cell_barcode_list = colnames(Seurat_obj_pass)
################################################

########################################################
########## HTO_DSB 
########################################################

############################### Read in RAW data for Empty droplet processing ############################
counts_FULL   <- Read10X_h5(filename= paste0(cellranger_Dir,"raw_peak_bc_matrix.h5"))

barcode_common_2 = Reduce(intersect, list(colnames(counts_FULL), rownames(HTO_data), rownames(ADT_data), rownames(metadata)))
counts_FULL = counts_FULL[, barcode_common_2]
metadata_F  = metadata[barcode_common_2,]
HTO_F       = HTO_data[barcode_common_2,]
ADT_F       = ADT_data[barcode_common_2,]

Seurat_obj_DSB <- CreateSeuratObject(counts = counts_FULL,assay = 'peaks',
                                     project = 'ATAC',min.cells = 1, meta.data = metadata_F)
Seurat_obj_DSB[["HTO"]] <- CreateAssayObject(counts = t(HTO_F))
Seurat_obj_DSB[["ADT"]] <- CreateAssayObject(counts = t(ADT_F))
Seurat_obj_DSB          <- NormalizeData(Seurat_obj_DSB, assay = "HTO", normalization.method = "CLR")
Seurat_obj_DSB          <- ScaleData(Seurat_obj_DSB, assay = "HTO")
Seurat_obj_DSB          <- HTODemux(Seurat_obj_DSB,  assay = "HTO", positive.quantile = 0.999)

table(Seurat_obj_DSB$HTO_classification.global)
### Doublet Negative  Singlet
###   18362    39741     9859

Seurat_obj_DSB$AFTER_QC = is.element(colnames(Seurat_obj_DSB),cell_barcode_list)

table(Seurat_obj_DSB$AFTER_QC,Seurat_obj_DSB$HTO_classification.global)
#        Doublet Negative Singlet
#  FALSE   11557    39741    9858
#  TRUE     6805        0       1

###########################################

# empty droplets from HTO demultiplex as idents = "Negative"
Idents(Seurat_obj) = "HTO_classification.global"
NegativeObject   = subset(Seurat_obj_DSB, idents = "Negative", subset = AFTER_QC =="FALSE")
SNG_SeuratObject = subset(Seurat_obj_DSB, subset = AFTER_QC=="TRUE")

# before normalizing make sure the negative droplets from hashing are true negatives. 
# many ways to do this. Here is a simple heuristic:
NegativeObject_2 = subset(NegativeObject, subset = nCount_peaks < 40) ###############

# omit outliers which have too much ADT read counts. ############################# therhold == 0.999
positive_ADT_tmp    = GetAssayData(SNG_SeuratObject, assay = "ADT", slot = 'counts') %>% as.matrix()
threshold_isotype   = quantile(positive_ADT_tmp["Mouse-IgG1",],.999) %>% round()
outlier_ADT         = colnames(positive_ADT_tmp)[positive_ADT_tmp["Mouse-IgG1",]>=threshold_isotype]
data.frame(ADT_outlier=outlier_ADT) %>% 
 write.table(.,paste0("../output/Signac/after_DSB/ADT_isotype_outlier.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)
SNG_SeuratObject_2  = SNG_SeuratObject[,!is.element(colnames(SNG_SeuratObject),outlier_ADT)]

###################################################
# get the empty droplet matrix and singlet matrix
neg_ADT_matrix      = GetAssayData(NegativeObject_2,   assay = "ADT", slot = 'counts') %>% as.matrix()
positive_ADT_matrix = GetAssayData(SNG_SeuratObject_2, assay = "ADT", slot = 'counts') %>% as.matrix()

# normalize the data with dsb normalization 
isotypes = "Mouse-IgG1"
ADT_DSB_data = DSBNormalizeProtein(cell_protein_matrix = positive_ADT_matrix,
                                   empty_drop_matrix = neg_ADT_matrix, 
                                   use.isotype.control = TRUE, 
                                   isotype.control.name.vec = isotypes) %>% t() %>% as.data.frame() 

ADT_DSB_data %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,paste0("../output/Signac/after_DSB/ADT_count_afterDSB_data.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)

###################################################
neg_HTO_matrix      = GetAssayData(NegativeObject_2, assay = "HTO", slot = 'counts') %>% as.matrix()
positive_HTO_matrix = GetAssayData(SNG_SeuratObject_2, assay = "HTO", slot = 'counts') %>% as.matrix()

# normalize the data with dsb normalization 
HTO_DSB_data = DSBNormalizeProtein(cell_protein_matrix = positive_HTO_matrix,
                                   empty_drop_matrix   = neg_HTO_matrix, 
                                   use.isotype.control = FALSE) %>% t() %>% as.data.frame()

HTO_DSB_data %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,paste0("../output/Signac/after_DSB/HTO_count_afterDSB_data.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)

#######################################################

## HTO_ADT after DSB data
HTO_DSB_data = read.table(paste0("../output/Signac/after_DSB/HTO_count_afterDSB_data.txt"),header=T,row.names=1,stringsAsFactors=F)
ADT_DSB_data = read.table(paste0("../output/Signac/after_DSB/ADT_count_afterDSB_data.txt"),header=T,row.names=1,stringsAsFactors=F)

##
HTO_DSB_data_2 = apply(HTO_DSB_data,c(1,2),round)
HTO_DSB_data_2[HTO_DSB_data_2<0] = 0

Seurat_obj_pass_2 = Seurat_obj_pass[,!is.element(colnames(Seurat_obj_pass),outlier_ADT)]

Seurat_obj_pass_2[["HTODSB"]] = CreateAssayObject(counts = t(HTO_DSB_data_2)[,colnames(Seurat_obj_pass_2)])
Seurat_obj_pass_2[["ADTDSB"]] = CreateAssayObject(counts = t(ADT_DSB_data)[,colnames(Seurat_obj_pass_2)])

Seurat_obj_pass_2 = NormalizeData(Seurat_obj_pass_2, assay = "HTODSB", normalization.method = "CLR")
Seurat_obj_pass_2 = ScaleData(Seurat_obj_pass_2, assay = "HTODSB")
Seurat_obj_pass_2 = HTODemux(Seurat_obj_pass_2, assay = "HTODSB", positive.quantile = 0.9)

###################
### Seurat_obj_pass_2 = NormalizeData(Seurat_obj_pass_2, assay = "ADTDSB", normalization.method = "CLR")
### Seurat_obj_pass_2 = ScaleData(Seurat_obj_pass_2, assay = "ADTDSB")

###################
# Global classification results
#### table(Seurat_obj_pass_2$HTODSB_classification.global) %>% print()
#### table(Seurat_obj_pass_2$HTODSB_maxID) %>% print()
#### table(Seurat_obj_pass_2$HTODSB_classification) %>% print()

###################
HTO_convert = data.frame(HTO_classification=Seurat_obj_pass_2$HTODSB_classification)
HTO_convert$HTO = HTO_convert$HTO_classification

## only combination of sgRNA[NTC/CD4/CD3E/ZAP70/NFKB2] & sgGuide(Guide1/Guide2) were picked up
## others are set as Negative
HTO_convert$HTO[grep("sgGuide",HTO_convert$HTO_classification,invert=T)] = "Negative"
HTO_convert$HTO[grep("sgGuide1",HTO_convert$HTO_classification)] = HTO_convert$HTO_classification[grep("sgGuide1",HTO_convert$HTO_classification)] %>%
                                                                                 gsub("_sgGuide1","",.) %>% gsub("sgGuide1_","",.) %>% paste0(.,"_sgGuide1")
HTO_convert$HTO[grep("sgGuide2",HTO_convert$HTO_classification)] = HTO_convert$HTO_classification[grep("sgGuide2",HTO_convert$HTO_classification)] %>%
                                                                                 gsub("_sgGuide2","",.) %>% gsub("sgGuide2_","",.) %>% paste0(.,"_sgGuide2")
## For this experiment, sgCD3E+sgGuide1 is sgCD3E CD4 DKO. Therefore, terminology are changed manually.
HTO_convert$HTO[HTO_convert$HTO_classification=="sgCD3E_sgGuide1"] = "sgCD3ECD4_sgGuide1"

### reading color list
sgRNA_list = read.table("data/sgRNA_list.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")
##                   HTO   color
## 1      sgNTC_sgGuide1 #E76BF3
## 2      sgNTC_sgGuide2 #E76BF3
## 3      sgCD4_sgGuide1 #00B0F6
## 4      sgCD4_sgGuide2 #00B0F6
## 5     sgCD3E_sgGuide2 #00BF7D
## 6  sgCD3ECD4_sgGuide1 #3ACCBD
## 7    sgZAP70_sgGuide1 #A3A500
## 8    sgZAP70_sgGuide2 #A3A500
## 9    sgNFKB2_sgGuide1 #F8766D
## 10   sgNFKB2_sgGuide2 #F8766D


Seurat_obj_pass_2$HTO_final = factor(HTO_convert$HTO,levels=c(sgRNA_list$HTO,"Negative"))
Seurat_obj_f = subset(x = Seurat_obj_pass_2,
                      subset = HTO_final!="Negative")
Seurat_obj_f$HTO_final           = factor(Seurat_obj_f$HTO_final,levels=sgRNA_list$HTO           )
Seurat_obj_f$HTO_final_Ridgeplot = factor(Seurat_obj_f$HTO_final,levels=sgRNA_list$HTO %>% rev(.))
## Seurat_obj_f$HTO_classification_final_2 = as.character(Seurat_obj_f$HTO_classification_final) %>% factor(.,levels=feature_list_2)

###################
# Global classification results
table(Seurat_obj_f$HTO_final) %>% print()
##    sgNTC_sgGuide1     sgNTC_sgGuide2     sgCD4_sgGuide1     sgCD4_sgGuide2
##               560                556                548                620
##   sgCD3E_sgGuide2 sgCD3ECD4_sgGuide1   sgZAP70_sgGuide1   sgZAP70_sgGuide2
##               621                639                614                639
##  sgNFKB2_sgGuide1   sgNFKB2_sgGuide2
##               534                494

print(Seurat_obj_f)
## An object of class Seurat
## 94611 features across 5825 samples within 5 assays
## Active assay: peaks (94523 features, 0 variable features)
##  4 other assays present: HTO, ADT, HTODSB, ADTDSB

########################################################
########## RidgePlot of HTO / ADT 
########################################################

# Group cells based on the max HTO signal
Idents(Seurat_obj_f) <- "HTO_final_Ridgeplot"
col_list = sgRNA_list$color %>% rev()

##
p1 = RidgePlot(Seurat_obj_f, assay = "HTO", features = HashTag_list, ncol = 7,slot = "data",same.y.lims = TRUE,cols=col_list)
pdf("../output/Signac/RidgePlot/HTO_RidgePlot.pdf",h=5,w=35)
  (p1) + plot_annotation( title = paste0(SAMPLE_TMP," after filter HTO"))
dev.off()

##
p1 = RidgePlot(Seurat_obj_f, assay = "ADT", features = rownames(Seurat_obj_f[["ADT"]]), ncol = 7,slot = "data",same.y.lims = TRUE,cols=col_list)
pdf("../output/Signac/RidgePlot/ADT_RidgePlot.pdf",h=30,w=35)
  (p1) + plot_annotation( title = paste0(SAMPLE_TMP," after filter ADT"))
dev.off()


########################################################
########## VlnPlot of HTO / ADT 
########################################################

Idents(Seurat_obj_f) <- "HTO_final"

##
p1 = VlnPlot(Seurat_obj_f, assay = "HTO", features = HashTag_list, ncol = 7,slot = "scale.data",same.y.lims = TRUE,cols=sgRNA_list$color, pt.size = 0)
pdf("../output/Signac/VlnPlot/HTO_VlnPlot.pdf",h=5,w=35)
  (p1) + plot_annotation( title = paste0(SAMPLE_TMP," after filter HTO"))
dev.off()

##
p1 = VlnPlot(Seurat_obj_f, assay = "ADT", features = rownames(Seurat_obj_f[["ADT"]]), ncol = 7,slot = "scale.data",same.y.lims = TRUE,cols=sgRNA_list$color, pt.size = 0)
pdf("../output/Signac/VlnPlot/ADT_VlnPlot.pdf",h=30,w=35)
  (p1) + plot_annotation( title = paste0(SAMPLE_TMP," after filter ADT"))
dev.off()

########################################################
### Saving Data After Filter                         ###
########################################################

##########
ADT_data          =  GetAssayData(Seurat_obj_f, assay = "ADT",   slot = "data") %>% t() %>% as.data.frame()
HTO_data          =  GetAssayData(Seurat_obj_f, assay = "HTO",   slot = "data") %>% t() %>% as.data.frame()
ADT_DSB_data      =  GetAssayData(Seurat_obj_f, assay = "ADTDSB",slot = "data") %>% t() %>% as.data.frame()
HTO_DSB_data      =  GetAssayData(Seurat_obj_f, assay = "HTODSB",slot = "data") %>% t() %>% as.data.frame()
HTO_res           =  data.frame(HTO = Seurat_obj_f$HTO_final)
HTO_res$HTO_sgRNA =  HTO_res$HTO %>% as.character() %>% strsplit(.,"_") %>% sapply(.,function(x){x[1]})

HTO_res         %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,"../output/Signac/after_filter_Signac/HTO_res_filtered.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)
HTO_data        %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,"../output/Signac/after_filter_Signac/HTO_filtered.txt"    ,row.names=F, col.names=T, sep="\t", append=F, quote=F)
ADT_data        %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,"../output/Signac/after_filter_Signac/ADT_filtered.txt"    ,row.names=F, col.names=T, sep="\t", append=F, quote=F)
HTO_DSB_data    %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,"../output/Signac/after_filter_Signac/HTO_DSB_filtered.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)
ADT_DSB_data    %>% rownames_to_column("CellBarcode") %>% 
 write.table(.,"../output/Signac/after_filter_Signac/ADT_DSB_filtered.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#**********************************************
## HTO heatmap _DSB soore
#**********************************************

HTO_DSB_2    = HTO_DSB_data
HTO_DSB_2[HTO_DSB_2<0.5] = 0    ######## too low expression is condsidered to be noise
HTO_DSB_2[HTO_DSB_2>2]   = 2    ########

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

data = HTO_DSB_2[c(hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[1] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[2] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[3] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[4] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[5] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[6] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[7] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[8] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[9] ],]),
                   hclust_order(HTO_DSB_2[rownames(HTO_res)[HTO_res$HTO==sgRNA_list$HTO[10] ],]) ),
                   c("sgNTC","sgCD4","sgCD3E","sgZAP70","sgNFKB2","sgGuide1","sgGuide2")] %>% 
       t() %>% as.data.frame()

col_list_HTOheatmap = colorRampPalette(brewer.pal(9,"Blues"))

df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)

df2 = left_join(data.frame(CellBarcode=colnames(data)),HTO_res%>%rownames_to_column("CellBarcode"),by="CellBarcode")
df2$HTO = factor(df2$HTO,levels=sgRNA_list$HTO)

p1  = ggplot(df,aes(x=colnames,y=rownames,fill=score))+
       ggrastr::geom_tile_rast()+
       scale_fill_gradientn(colours=col_list_HTOheatmap(250))+ 
       theme_bw() +
       scale_x_discrete(limits=unique(df$colnames)) +
       scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
       theme(axis.text.x= element_blank(),
             axis.text.y= element_text(colour="black",size=10,face="bold", family = "Helvetica"),
             axis.title = element_text(colour="black",size=10,face="bold", family = "Helvetica"),
             strip.text = element_text(size=8),
             axis.line  = element_line(colour = "black"),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
       labs(x="Cell Barcode",y="HTO")

p2 = ggplot(df2, aes(x=CellBarcode,y=1,fill=HTO))+
      ggrastr::geom_tile_rast()+
      scale_fill_manual(values = sgRNA_list$color)+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_discrete(limits=colnames(data)) +
      theme(axis.ticks  = element_blank(),
            axis.text   = element_blank(),
            axis.title.x= element_blank(),
            axis.title.y= element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=10,face="bold", family = "Helvetica"),
            plot.title  = element_text(size=10,face="bold", family = "Helvetica"),
            plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"),legend.position="none")+
      labs(title=paste0(SAMPLE_TMP," HTO afterDSB"),y="sgRNA")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)

maxWidth = grid::unit.pmax(gp1$widths, gp2$widths)
gp1$widths <- gp2$widths <- as.list(maxWidth)

################ Fig B ################################################################################
pdf("../plots/HTO_afterDSB_heatmap.pdf",h=5,w=10)
 grid.arrange(gp2,gp1, ncol=1,heights=c(0.1,0.9))
dev.off()

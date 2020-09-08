options(stringsAsFactors=F)

#******************************************
### args        = commandArgs(trailingOnly = T)
### SAMPLE_TMP  = as.character( args[1] )
SAMPLE_TMP  = "Perturb_CD4_stim"

##################
fragment_path = paste0("../../../asap_large_data_files/CD4_stimulation_asapseq/cellranger/fragments.tsv.gz")

suppressMessages(library("ArchR"))
suppressMessages(library("patchwork"))
suppressMessages(library("tidyverse"))
suppressMessages(library("gridExtra"))
suppressMessages(library("data.table"))

set.seed(111)

Single_cell_list     = read.table("../output/Signac/after_filter_Signac/HTO_res_filtered.txt",header=T,row.names=NULL,stringsAsFactors=F)
Single_cell_list$UMI = Single_cell_list$CellBarcode %>% paste0(SAMPLE_TMP,"#",.)

sgRNA_list     = read.table("../data/sgRNA_list.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")
genotype_info  = read.table("../output/genotype/final_res/CellBarcode_genotype_result.txt",header=T,row.names=NULL,stringsAsFactors=F,sep="\t")
colnames(genotype_info) = c("CellBarcode","genotype")
genotype_info$genotype = paste0("Genotype_",genotype_info$genotype)

#dir_list = c("output/ArchR/arrow/","output/ArchR/ArchRLogs","output/ArchR/QC","output/ArchR/chromVAR/WMW_res","output/ArchR/DimPlot","output/ArchR/DiffPeaks","output/ArchR/Rdata")
#for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}

#****************************************#
# Making of ArchR project
#****************************************#

addArchRThreads(threads = 30)
addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles      = fragment_path ,
  sampleNames     = SAMPLE_TMP ,
  filterTSS       = 1, #Dont set this too high because you can always increase later
  filterFrags     = 500, 
  addTileMat      = TRUE,
  addGeneScoreMat = FALSE,
  force = TRUE,
  QCDir = "../output/ArchR/QC",
  logFile = createLogFile(name = "createArrows", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging())
)

doubScores <- addDoubletScores(input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force = TRUE,
  outDir = "../output/ArchR/QC",
  logFile = createLogFile(name = "addDoubletScores", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging())
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "../output/ArchR",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

paste0(SAMPLE_TMP,".arrow") %>% file.remove(.)

## cell barcode after QC filter by Seurat
target_cell = intersect(Single_cell_list$UMI,rownames(proj))
proj_subset = subsetCells(ArchRProj = proj, cellNames = target_cell)

Single_cell_list = left_join(data.frame(UMI=target_cell),Single_cell_list,by="UMI")
proj_subset$HTO       = Single_cell_list$HTO
proj_subset$HTO_sgRNA = Single_cell_list$HTO_sgRNA

fragment_info  = getCellColData(proj_subset, select = "nFrags") %>% as.data.frame(.)
fragment_info %>% rownames_to_column("UMI") %>%
 write.table(.,"../output/ArchR/QC/fragment_info.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)
summary(fragment_info) %>% print(.)
##      nFrags
##  Min.   : 3535
##  1st Qu.: 9599
##  Median :14731
##  Mean   :17332
##  3rd Qu.:22426
##  Max.   :59892

#****************************************#
# Dimension Reduction by UMAP
#****************************************#

proj_subset <- addIterativeLSI(
    ArchRProj = proj_subset,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2,0.8), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 12000, 
    dimsToUse = 1:20,
    force=TRUE,
    logFile = createLogFile(name = "addIterativeLSI", logDir = "../output/ArchR/ArchRLogs", useLogs = getArchRLogging())
)

proj_subset <- ArchR::addClusters(
    input = proj_subset,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force=TRUE,
    logFile = createLogFile(name = "addClusters", logDir = "../output/ArchR/ArchRLogs", useLogs = getArchRLogging())
)

proj_subset <- addUMAP(
    ArchRProj = proj_subset, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 20, 
    minDist = 0.5, 
    metric = "cosine",
    force=TRUE
)

#****************************************#
# Drawing DimPlot using UMAP info
#****************************************#

#### plot figures using ggplot2 (all sgRNA information)
UMAP_res     = proj_subset@embeddings[["UMAP"]][[1]]
colnames(UMAP_res) = c("UMAP_1","UMAP_2")
UMAP_res %>% rownames_to_column("UMI") %>% 
 write.table(.,"../output/ArchR/DimPlot/ArchR_UMAP_res.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

UMAP_res_2 = UMAP_res %>% rownames_to_column("UMI") %>% 
              left_join(.,Single_cell_list,by="UMI")

UMAP_res_2$HTO  = factor(UMAP_res_2$HTO,levels=sgRNA_list$HTO)

p1 = ggplot(UMAP_res_2, aes(x=UMAP_1,y=UMAP_2,color=HTO)) + 
      geom_point(size = .5) + 
      theme_void() + 
      scale_colour_manual(values=sgRNA_list$color)+
      labs(x="UMAP_1",y="UMAP_2")
p2 = AugmentPlot(p1, h=5, w=5, dpi=100) + plot_annotation(title = paste0(SAMPLE_TMP," varFeature 12000 PC 1:20 DimPlot with ATACpeak"))

paste0("../output/ArchR/DimPlot/ArchR_UMAP_varFeature12000_PC_1_20_DimPlot_v2.pdf") %>% pdf(.,h=5,w=5)
 p2
dev.off()

#### plot figures using ggplot2 (each sgRNA information)
dim_Plot = list()
for(tmp_num in 1:nrow(sgRNA_list)){
    dim_Plot[[tmp_num]] = ggplot(UMAP_res_2, aes(x=UMAP_1,y=UMAP_2))+
            geom_point(size = .5,col="gray85") +
            geom_point(data = UMAP_res_2 %>% dplyr::filter(HTO==sgRNA_list$HTO[tmp_num]),aes(x=UMAP_1,y=UMAP_2),size = 1,col=sgRNA_list$color[tmp_num]) +
            theme_void() +
            labs(x="UMAP_1",y="UMAP_2")
    dim_Plot[[tmp_num]] = AugmentPlot(dim_Plot[[tmp_num]], h=5, w=5, dpi=100) 
}

### Fig B #################################################################################################
paste0("../plots/ArchR_UMAP_DimPlot.pdf") %>% pdf(.,h=10,w=25)
 (dim_Plot[[1]]|dim_Plot[[3]]|dim_Plot[[5]]|dim_Plot[[7]]|dim_Plot[[9]])/(dim_Plot[[2]]|dim_Plot[[4]]|dim_Plot[[6]]|dim_Plot[[8]]|dim_Plot[[10]]) + 
  plot_annotation(title = paste0(SAMPLE_TMP," varFeature 12000 PC 1:20 DimPlot with ATACpeak"))
dev.off()

#****************************************#
# Detecting Diffeential Peaks
#****************************************#

proj_subset = addGroupCoverages(ArchRProj = proj_subset, groupBy = "HTO",
                                logFile   = createLogFile(name = "addGroupCoverages", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))
pathToMacs2 = findMacs2()
proj_subset = addReproduciblePeakSet(ArchRProj = proj_subset,
                                     groupBy = "HTO",
                                     pathToMacs2 = pathToMacs2,
                                     logFile = createLogFile(name = "addReproduciblePeakSet", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))

proj_subset = addPeakMatrix(proj_subset,
                            logFile   = createLogFile(name = "addPeakMatrix", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))

markersPeaks <- getMarkerFeatures(ArchRProj = proj_subset, 
                                  useMatrix = "PeakMatrix", 
                                  groupBy = "HTO",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  testMethod = "wilcoxon",
                                  logFile = createLogFile(name = "getMarkerFeatures", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()) )
markerList   <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5")
heatmapPeaks <- markerHeatmap(seMarker = markersPeaks, 
                              cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5",
                              transpose = TRUE,
                              logFile = createLogFile(name = "plotMarkerHeatmap", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()) )

pdf("../output/ArchR/DiffPeaks/heatmapPeaks_markerPeak.pdf",h=5,w=5)
 draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()

bed_peak      = rowData(markersPeaks) %>% as.data.frame() %>% dplyr::select(seqnames,start,end)
bed_peak$name = formatC(1:nrow(bed_peak),width=6,flag="0") %>% paste0("Peak_",.)
write.table(bed_peak,"../output/ArchR/DiffPeaks/Union_peak.bed",row.names=F, col.names=F, sep="\t", append=F, quote=F)

#####################

#****************************************#
# chromVAR  ## background is set to NTC 
#****************************************#

proj_subset  = addMotifAnnotations(ArchRProj = proj_subset , motifSet = "cisbp", name = "Motif",
                                   logFile = createLogFile(name = "addMotifAnnotations", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))
proj_subset  = addArchRAnnotations(ArchRProj = proj_subset , collection = "EncodeTFBS",
                                   logFile = createLogFile(name = "addArchRAnnotations", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))

# get peak x cell matrix
counts_SE  <- getMatrixFromProject(proj_subset,"PeakMatrix",
                                   logFile = createLogFile(name = "getMatrixFromProject", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))
counts_mat <- assays(counts_SE)[["PeakMatrix"]]
counts_mat <- counts_mat[,Single_cell_list$UMI[Single_cell_list$HTO_sgRNA=="sgNTC"]]

##--- MODIFY HERE-- ##
# subset the rows of counts_mat to only NTC cells
rS = data.frame(proj_subset@peakSet,
                rowSums = Matrix::rowSums(counts_mat) )

se <- SummarizedExperiment::SummarizedExperiment(
  assays = SimpleList(counts = as.matrix(data.frame(rS$rowSums, 1))),
  rowData = DataFrame(bias = rS$GC)
)

bgdPeaks <- chromVAR::getBackgroundPeaks(
  object = se,
  bias = rowData(se)$bias
)

bgdPeaks = SummarizedExperiment(assays = SimpleList(bgdPeaks = bgdPeaks), 
                                 rowRanges = GRanges(rS$seqnames,IRanges(rS$start,rS$end),value=rS$rowSums,GC=rS$GC))
proj_subset = addDeviationsMatrix(bgdPeaks = bgdPeaks,  
                                  ArchRProj = proj_subset, 
                                  peakAnnotation = "Motif",
                                  force = TRUE,
                                  logFile = createLogFile(name = "addDeviationsMatrix", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))

plotVarDev = getVarDeviations(proj_subset, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj_subset, addDOC = FALSE)

MotifMatrix_tmp = getMatrixFromProject(ArchRProj = proj_subset,
                                       useMatrix = "MotifMatrix",
                                       useSeqnames = NULL,
                                       verbose = TRUE,
                                       binarize = FALSE,
                                       threads = 100,
                                       logFile = createLogFile(name = "getMatrixFromProject", logDir = "output/ArchR/ArchRLogs", useLogs = getArchRLogging()))
zscore_matrix     = assays(MotifMatrix_tmp)$z          %>% as.data.frame()
deviations_matrix = assays(MotifMatrix_tmp)$deviations %>% as.data.frame()

zscore_matrix     %>% rownames_to_column("motif") %>%
 write.table(.,"../output/ArchR/chromVAR/motif_zscore_ChromVAR.txt"    ,row.names=F, col.names=T, sep="\t", append=F, quote=F)
deviations_matrix %>% rownames_to_column("motif") %>%
 write.table(.,"../output/ArchR/chromVAR/motif_deviations_ChromVAR.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#****************************************#
# Drawing Heatmap using chromVAR info
#****************************************#

#### Results of several motif enrichment look almost same, therefore some motifs are omitted manually.
omit_TF = c("ENSG00000232923_590","ENSG00000233294_591","ENSG00000233361_592","ENSG00000233698_593","ENSG00000233794_594","ENSG00000234550_596","ENSG00000236207_598","ENSG00000237516_601","ENSG00000230776_587","ENSG00000215754_570","ENSG00000224558_572","ENSG00000225958_574")

VarDeviation = getVarDeviations(proj_subset, name = "MotifMatrix", plot = FALSE) %>% as.data.frame()
write.table(VarDeviation,paste0("../output/ArchR/chromVAR/chromVAR_deviation_rank.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)

#### Visualization of Top 100 TFs
TF_rank_num = 100
Top_TF_target  = setdiff(VarDeviation$name,omit_TF) %>% head(.,n=TF_rank_num)

zscore_part = zscore_matrix[Top_TF_target,]
zscore_part[zscore_part>2]  = 2
zscore_part[-2>zscore_part] = -2

colname_list = left_join(data.frame(UMI=colnames(zscore_part)),Single_cell_list, by = "UMI" ) %>%
                left_join(. , genotype_info, by = "CellBarcode")
colname_list$HTO   = factor(colname_list$HTO, levels=sgRNA_list$HTO)
colname_list       = colname_list %>% arrange(HTO,genotype)
colname_list$group = paste0(colname_list$HTO,"_genotype_",colname_list$genotype)
group_list         = unique(colname_list$group)

#### clustring was done in each sgRNA. (with ward.D2)
hclust_order = function(DATA,METHOD="ward.D2"){
  suppressMessages(library(ggdendro))
  rd     = dist(DATA)
  hc     = hclust(d=rd,method=METHOD)
  dhc    = as.dendrogram(hc)
  ddata  = dendro_data(dhc, type = "rectangle")
  col_order = as.character(ddata$labels$label)
  col_order
}

#### reoderding in each sgRNA/genotype
colname_order_list = list()
for(i in 1:length(group_list)){
  tmp_data = zscore_part[,colname_list$UMI[colname_list$group==group_list[i]]]
  colname_order_list[[i]] = tmp_data %>% t(.) %>% hclust_order(.)
}

colname_list_f       = left_join(data.frame(UMI=unlist(colname_order_list)),Single_cell_list, by = "UMI" ) %>%
                        left_join(. , genotype_info, by = "CellBarcode")
colname_list_f$HTO   = factor(colname_list_f$HTO, levels=sgRNA_list$HTO)

data = zscore_part[hclust_order(zscore_part),colname_list_f$UMI]

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
col_list   = colorRampPalette(blueYellow)

df = gather(data.frame(rownames=rownames(data),data),key=colnames,value=score,-rownames)

p1  = ggplot(df,aes(x=colnames,y=rownames,fill=score))+
    ggrastr::geom_tile_rast()+
    scale_fill_gradientn(colours=col_list(250))+
    theme_bw() +
    scale_x_discrete(limits=unique(df$colnames)) +
    scale_y_discrete(limits=unique(df$rownames)%>%rev()) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(colour="black",size=1.5,face="bold", family = "Helvetica"),
          axis.title=element_text(colour="black",size=5,face="bold", family = "Helvetica"),
          strip.text=element_text(size=8),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(x="cell_barcode",y="TF (top 100)")

## for(i in 1:9){p1 = p1 + geom_vline(xintercept=tmp_sgRNA_num[i]+0.5,col="black") }

p2 = ggplot(colname_list_f, aes(x=UMI,y=1,fill=HTO))+
      ggrastr::geom_tile_rast()+scale_fill_manual(values = sgRNA_list$color )+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_discrete(limits=colname_list_f$UMI) +
      theme(axis.title.x=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
            axis.text.x=element_blank(),axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=5,face="bold", family = "Helvetica"),
            plot.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
            plot.margin = unit(c(0.2,0.1,0.2,0.1), "cm"),legend.position="none")+
      labs(y="sgRNA")

p3 = ggplot(colname_list_f, aes(x=UMI,y=1,fill=genotype))+
      ggrastr::geom_tile_rast()+scale_fill_manual(values = c("#D58DB6","#90C37A","#F5D18E") )+
      scale_y_continuous(expand=c(0,0)) +
      scale_x_discrete(limits=colname_list_f$UMI) +
      theme(axis.title.x=element_blank(),axis.ticks=element_blank(),axis.text.y=element_blank(),
            axis.text.x=element_blank(),axis.title.y=element_text(colour="black",angle = 0,vjust=0.5,hjust=1,size=5,face="bold", family = "Helvetica"),
            plot.title=element_text(colour="black",size=10,face="bold", family = "Helvetica"),
            plot.margin = unit(c(0.2,0.1,0.2,0.1), "cm"),legend.position="none")+
      labs(y="genotype")

gp1 = ggplotGrob(p1)
gp2 = ggplotGrob(p2)
gp3 = ggplotGrob(p3)

maxWidth = grid::unit.pmax(gp1$widths, gp2$widths, gp3$widths)
gp1$widths <- gp2$widths <- gp3$widths <- as.list(maxWidth)

### Fig D #################################################################################################
paste0("../plots/chromVAR_heatmap_Top100_motif.pdf") %>% pdf(.,h=3,w=5)
 grid.arrange(gp2,gp3,gp1, ncol=1,heights=c(0.1,0.1,0.9))
dev.off()

data %>% rownames_to_column("motif") %>%
 write.table(.,paste0("../output/ArchR/chromVAR/chromVAR_top100_for_plot.txt"),row.names=F, col.names=T, sep="\t", append=F, quote=F)

#****************************************#
# Drawing VolcanoPlot using chromVAR info -> Wilcoxon-Mann-Whitney
#****************************************#

deviations_matrix = fread("../output/ArchR/chromVAR/motif_deviations_ChromVAR.txt") %>% as.data.frame(.) %>% column_to_rownames("motif")
zscore_matrix     = fread("../output/ArchR/chromVAR/motif_zscore_ChromVAR.txt") %>% as.data.frame(.) %>% column_to_rownames("motif")

#### Extracting deviation score
for(sgRNA_tmp in sgRNA_list$HTO[3:10]){
    barcode_tmp_1 = Single_cell_list$UMI[is.element(Single_cell_list$HTO,c(sgRNA_list$HTO[1:2]))]
    barcode_tmp_2 = Single_cell_list$UMI[is.element(Single_cell_list$HTO,c(sgRNA_tmp))]
    df_tmp_1      = deviations_matrix[,barcode_tmp_1]
    df_tmp_2      = deviations_matrix[,barcode_tmp_2]

    for(i in 1:nrow(deviations_matrix)){
      TF_name  = rownames(deviations_matrix)[i]
      df_tmp_3 = data.frame(sgRNA=c(rep("sgNTC",ncol(df_tmp_1)),rep(sgRNA_tmp,ncol(df_tmp_2))),
                            deviation=c(df_tmp_1[i,] %>% as.numeric(.),df_tmp_2[i,] %>% as.numeric(.)) )
      dir.create(paste0("output/ArchR/chromVAR/WMW_raw_data/",sgRNA_tmp), showWarnings = FALSE, recursive = TRUE)
      paste0("output/ArchR/chromVAR/WMW_raw_data/",sgRNA_tmp,"/",SAMPLE_TMP,"_",sgRNA_tmp,"_",TF_name,"_raw_data.txt") %>% 
       write.table(df_tmp_3,.,row.names=F, col.names=T, sep="\t", append=F, quote=F)
    }
}

#### Extracting median of deviation Zscore
df_tmp      = zscore_matrix[,Single_cell_list$UMI[is.element(Single_cell_list$HTO,c(sgRNA_list$HTO[1:2]))]] %>% 
               apply(.,1,median)
res_tmp     = data.frame(TF = names(df_tmp),sgRNA=df_tmp)
colnames(res_tmp)[2] = "sgNTC"
res_sum = res_tmp

for(sgRNA_tmp in sgRNA_list$HTO[3:10]){
    df_tmp      = zscore_matrix[,Single_cell_list$UMI[is.element(Single_cell_list$HTO,sgRNA_tmp)]] %>% 
                   apply(.,1,median)
    res_tmp     = data.frame(TF = names(df_tmp),sgRNA=df_tmp)
    colnames(res_tmp)[2] = sgRNA_tmp
    res_sum = full_join(res_sum,res_tmp,by="TF")
}

write.table(res_sum,"../output/ArchR/chromVAR/WMW_res/deviation_zScore_median_data.txt",row.names=F, col.names=T, sep="\t", append=F, quote=F)

#****************************************#
# Saving information
#****************************************#

save(list=ls(),file=paste0("../../../asap_large_data_files/CD4_stimulation_asapseq/output/ArchR/ArchR.Rdata"))

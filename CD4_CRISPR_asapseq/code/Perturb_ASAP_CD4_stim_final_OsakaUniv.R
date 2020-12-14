#-------------------------------------------------------------------#
# Perturb CD4 stimulation log final _ for Github
# 200813 written by Y.Takeshima / Kelvin Chen @ Osaka.Univ, JPN
#-------------------------------------------------------------------#

### setting of cutrrent working directory
###  WD=Perturb_ASAP_CD4_stim_final/asap_reproducibility/CD4_stimulation_asapseq; ##############
###  cd ${WD}; ##############

WD=$(pwd);
SAMPLE_TMP=Perturb_CD4_stim;

#****************************************#
# Seurat Signac 
#****************************************#

JOB_NAME=${SAMPLE_TMP}_Signac;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R0_Signac.sh \
      code/11_ATACseq_Signac.R \
      ${SAMPLE_TMP} \
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

#****************************************#
# genotyping
# ref) https://github.com/statgen/popscle.git
# ref) https://github.com/aertslab/popscle_helper_tools
#****************************************#

JOB_NAME=${SAMPLE_TMP}_freebayes_popscle_freemuxlet;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh code/21_freebayes_popscle_freemuxlet.sh \
 ${SAMPLE_TMP}\
 1> log/${JOB_NAME}/freebayes_popscle_freemuxlet.log \
 2> error/${JOB_NAME}/freebayes_popscle_freemuxlet.error &

#****************************************#
# ArchR
#****************************************#

JOB_NAME=${SAMPLE_TMP}_ArchR_1;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/31_ArchR.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

## Wilcoxon Mann Whitnney analysis
for sgRNA_TMP in $(ls output/ArchR/chromVAR/WMW_raw_data);do
 JOB_NAME=${SAMPLE_TMP}_WilcoxonMannWhitnney;
 mkdir -p log/${JOB_NAME} error/${JOB_NAME};
 DD_TMP=${WD}/output/ArchR/chromVAR/WMW_raw_data
 nohup sh code/32_WilcoxonMannWhitnney.sh \
       ${SAMPLE_TMP}\
       ${DD_TMP} \
       ${sgRNA_TMP} \
       1> log/${JOB_NAME}/${SAMPLE_TMP}_${sgRNA_TMP}.log \
       2> error/${JOB_NAME}/${SAMPLE_TMP}_${sgRNA_TMP}.error &
done

rm -r output/ArchR/chromVAR/WMW_raw_data

#****************************************#
# chromVAR volcanoPlot scatterPlot
#****************************************#
JOB_NAME=${SAMPLE_TMP}_chromVAR_volcano_scatter_plot;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/34_chromVAR_volcano_scatter_plot.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

#****************************************#
# chromVar dimPlot
#****************************************#
JOB_NAME=${SAMPLE_TMP}_chromVAR_Dimplot;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/35_chromVAR_Dimplot.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &


#****************************************#
# GeneScore --manual calculation--
#****************************************#

## divide bam file into each CellBarcode
## pick up target region BED.
mkdir -p output/GeneScore/job_list

awk 'OFS="\t"''{print $1,$2"_"$1}' output/genotype/cell_barcode_list.txt \
| split -d -a 3 -l 100 - output/GeneScore/job_list/cell_barcode_list_split;

for JOB_ID in $(seq -f %03g 0 58);do  # 41 50 # 51 58       #0 58
 JOB_NAME=${SAMPLE_TMP}_sinto_CellBarcode;
 mkdir -p log/${JOB_NAME} error/${JOB_NAME};
 nohup sh code/41_sinto_CellBarcode.sh \
       ${JOB_ID} \
       1> log/${JOB_NAME}/${SAMPLE_TMP}_split${JOB_ID}.log \
       2> error/${JOB_NAME}/${SAMPLE_TMP}_split${JOB_ID}.error &
done

rm -r output/GeneScore/bam  output/GeneScore/bedpe output/GeneScore/tmp_bam  output/GeneScore/tmp_bed  output/GeneScore/tmp_bedpe

## calculate each bed number fgramnents w/o Blacklist.
wc -l output/GeneScore/bed/*/*.bed \
| sed "s/^[ \t]*//" - \
| grep -v total - > output/GeneScore/bed_fragment_number_list.txt

## summarizing GeneScore 
## drawing scatterplot of GeneScore and ADT zScore
JOB_NAME=${SAMPLE_TMP}_GeneScore_plot;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/42_GeneScore_plot.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

rm -r output/GeneScore/bed output/GeneScore/bed_target output/GeneScore/job_list

#****************************************#
# ADT heatmap
#****************************************#
JOB_NAME=${SAMPLE_TMP}_ADT_hetamap;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/51_ADT_hetamap.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

#****************************************#
# ADT dimplot
#****************************************#
JOB_NAME=${SAMPLE_TMP}_ADT_Dimplot;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/52_ADT_Dimplot.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

#****************************************#
# each sgRNA bam -> MACS2 -> bigwig
#****************************************#

JOB_NAME=${SAMPLE_TMP}_sinto;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh code/61_sinto.sh \
      ${SAMPLE_TMP} \
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

######### Making of bigwig file in each sgRNA
for BAM_TMP in $(ls output/each_sgRNA/bam/*);do
    echo ${BAM_TMP}
    COND_TMP=$(basename ${BAM_TMP} | cut -d "." -f 1)
    JOB_NAME=${SAMPLE_TMP}_bam_sort_MACS2
    mkdir -p log/${JOB_NAME} error/${JOB_NAME}
    nohup sh code/62_bam_sort_MACS2_hg38.sh \
          ${BAM_TMP}\
          1> log/${JOB_NAME}/${COND_TMP}.log \
          2> error/${JOB_NAME}/${COND_TMP}.error &
done

rm -r output/each_sgRNA/bam output/each_sgRNA/MACS2

#****************************************#
# merge filtered peaks (w/o black list) in each sgRNA
#****************************************#
cat output/each_sgRNA/MACS2_peak/* \
| sort -k 1,1 -k 2n,2 - \
| bedtools merge -i - > output/each_sgRNA/MACS2_peak_call_merge.bed

#****************************************#
# motif match NFKB2
#****************************************#
JOB_NAME=${SAMPLE_TMP}_motifmatchr_NFKB;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh ~/shell/R1_Signac.sh \
      code/71_motifmatchr_NFKB.R \
      ${SAMPLE_TMP}\
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &


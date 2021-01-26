#---------------------------------------#
# Perturb24 CD4 stim Fig7/ FigS7
# 21015 written by Y.Takeshima / Kelvin Chen @ Osaka.Univ, JPN
#---------------------------------------#

### setting of cutrrent working directory
### cd ${WD};
SAMPLE_TMP=Perturb24_CD4_stim

#****************************************#
# Seurat -> ArchR
#****************************************

JOB_NAME=${SAMPLE_TMP}_Seurat_ASAP_UMAP;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh code/R1.sh  \
      code/11_Seurat_ArchR_Perturb24_CD4_stim.R \
      ${SAMPLE_TMP} \
      1> log/${JOB_NAME}/${JOB_NAME}.log \
      2> error/${JOB_NAME}/${JOB_NAME}.error &

#****************************************#
# Making of bigwig file
#****************************************

######### divide bam file sinto each sgRNA
JOB_NAME=${SAMPLE_TMP}_sinto;
mkdir -p log/${JOB_NAME} error/${JOB_NAME};
nohup sh code/21_sinto_Perturb24_CD4_stim.sh \
      ${SAMPLE_TMP} \
      1> log/${JOB_NAME}/${SAMPLE_TMP}.log \
      2> error/${JOB_NAME}/${SAMPLE_TMP}.error &

######### Making of bigwig file in each sgRNA
for BAM_TMP in $(ls output/sinto/bam_1/*.bam);do
    echo ${BAM_TMP}
    COND_TMP=$(basename ${BAM_TMP} | cut -d "." -f 1)
    JOB_NAME=${SAMPLE_TMP}_bam_sort_MACS2
    mkdir -p log/${JOB_NAME} error/${JOB_NAME}
    nohup sh code/22_bam_sort_MACS2_hg38_Perturb24_CD4_stim.sh \
          ${COND_TMP} \
          1> log/${JOB_NAME}/${COND_TMP}.log \
          2> error/${JOB_NAME}/${COND_TMP}.error &
done

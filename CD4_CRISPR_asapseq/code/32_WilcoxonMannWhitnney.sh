#!/bin/bash
#$ -S /bin/bash

SAMPLE_TMP=$1
DD_TMP=$2
sgRNA_TMP=$3

echo SAMPLE : ${SAMPLE_TMP}
echo DD : ${DD_TMP}
echo sgRNA : ${sgRNA_TMP}

######### 
for TF_TMP in $(ls output/ArchR/chromVAR/WMW_raw_data/${sgRNA_TMP} | cut -d "_" -f 6,7);do
    python \
     code/33_WMW.py \
     ${SAMPLE_TMP} \
     ${DD_TMP} \
     ${sgRNA_TMP} \
     ${TF_TMP} >> output/ArchR/chromVAR/WMW_res/${sgRNA_TMP}_WMW_summary.txt
done

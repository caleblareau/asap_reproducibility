#!/bin/bash
#$ -S /bin/bash

SAMPLE_TMP=$1
WD=$(pwd)

echo SAMPLE : ${SAMPLE_TMP}

######### Virtual environment
conda init bash            ##
source ~/.bash_profile     ##
conda activate PerturbASAP ##

################# making list
mkdir -p output/sinto/list output/sinto/bam_1 output/sinto/bam_2

grep "\-1" output/Signac/after_filter_Signac/Seurat_cluster_data.txt > output/sinto/list/HTO_convert_list_1.txt

grep "\-2" output/Signac/after_filter_Signac/Seurat_cluster_data.txt \
| sed -e 's/-2/-1/g' > output/sinto/list/HTO_convert_list_2.txt

######### dividing bam file into each sgRNA
for TMP_NUM in 1 2 ;do
 cd ${WD}/output/sinto/bam_${TMP_NUM} ;
 sinto filterbarcodes  -b ${WD}/data/cellranger/possorted_bam_${TMP_NUM}.bam \
                       -c ${WD}/output/sinto/list/HTO_convert_list_${TMP_NUM}.txt \
                       -p 30 ;
 cd ${WD}
done




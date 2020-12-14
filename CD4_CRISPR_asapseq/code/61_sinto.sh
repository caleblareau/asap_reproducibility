#!/bin/bash
#$ -S /bin/bash

SAMPLE_TMP=$1
WD=$(pwd)

echo SAMPLE : ${SAMPLE_TMP}

mkdir -p output/each_sgRNA/bam

######### dividing bam file into each sgRNA
cd ${WD}/output/each_sgRNA/bam ;
sinto filterbarcodes  -b ${WD}/../../asap_large_data_files/CD4_stimulation_asapseq/cellranger/possorted_bam.bam  \
                      -c ${WD}/output/genotype/cell_barcode_list.txt \
                      -p 30 ;
cd ${WD}

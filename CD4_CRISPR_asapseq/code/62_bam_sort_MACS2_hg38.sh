#!/bin/sh
#$ -S /bin/sh

########################################################
FILE_NAME=$1;
COND_TMP=$(basename ${FILE_NAME} | cut -d "." -f 1);

PEAK=output/each_sgRNA/MACS2/${COND_TMP}_peaks.narrowPeak;
FILTERED_PEAK=output/each_sgRNA/MACS2_peak/${COND_TMP}_peaks.narrowPeak.filt.bed

GENOME_FILE=data/UCSC_hg38_chromInfo.txt;
BLACKLIST=data/hg38-blacklist.v2.bed

echo COND_TMP : ${COND_TMP};

mkdir -p output/each_sgRNA/MACS2 output/each_sgRNA/bigwig output/each_sgRNA/MACS2_peak;

########################################################
samtools index ${FILE_NAME}

########################################################

##MACS2
smooth_window=150;
shiftsize=$((-$smooth_window/2));

macs2 callpeak \
 -t ${FILE_NAME}  \
 -f BAM -g hs \
 --outdir output/each_sgRNA/MACS2 \
 -n ${COND_TMP} \
 -B --nomodel \
 --extsize ${smooth_window} \
 --shift ${shiftsize} \
 --SPMR --call-summits --keep-dup 1\
 -q 0.1

bedGraphToBigWig output/each_sgRNA/MACS2/${COND_TMP}_treat_pileup.bdg ${GENOME_FILE} output/each_sgRNA/bigwig/${COND_TMP}_treat_pileup.bw

########## Blacklist filtering
bedtools intersect -v -a ${PEAK} -b ${BLACKLIST} \
| awk 'BEGIN {OFS="\t"} {if($5>1000)$5=1000;print $0}' \
| grep -P 'chr[0-9XY]+(?!_)'  > ${FILTERED_PEAK}

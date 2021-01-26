#!/bin/sh
#$ -S /bin/sh

########################################################
COND_TMP=$1

echo COND_TMP : ${COND_TMP}

########################################################
mkdir -p output/sinto/bam output/sinto/MACS2 output/sinto/bigwig

## hg38 Blacklist
GENOME_FILE=data/chromInfo.txt;
BLACKLIST=data/hg38-blacklist.v2.bed

## input data
PEAK=output/sinto/MACS2/${COND_TMP}_peaks.narrowPeak;
FILTERED_PEAK=output/sinto/MACS2/${COND_TMP}_peaks.narrowPeak.filt.bed

########################################################
samtools merge -@ 10 output/sinto/bam/${COND_TMP}.bam output/sinto/bam_1/${COND_TMP}.bam output/sinto/bam_2/${COND_TMP}.bam
samtools sort -@ 10 output/sinto/bam/${COND_TMP}.bam -o output/sinto/bam/${COND_TMP}_sort.bam
samtools index output/sinto/bam/${COND_TMP}_sort.bam

########################################################

##MACS2
smooth_window=150;
shiftsize=$((-$smooth_window/2))

~/.local/bin/macs2 callpeak \
 -t output/sinto/bam/${COND_TMP}_sort.bam \
 -f BAM -g hs \
 --outdir output/sinto/MACS2 \
 -n ${COND_TMP} \
 -B --nomodel \
 --extsize ${smooth_window} \
 --shift ${shiftsize} \
 --SPMR --call-summits --keep-dup 1\
 -q 0.1

~/tool/UCSC/bedGraphToBigWig output/sinto/MACS2/${COND_TMP}_treat_pileup.bdg ${GENOME_FILE} output/sinto/bigwig/${COND_TMP}_treat_pileup.bw

########## Blacklist filtering
/home/ytakeshima/.local/bin/bedtools intersect -v -a ${PEAK} -b ${BLACKLIST} \
| awk 'BEGIN {OFS="\t"} {if($5>1000)$5=1000;print $0}' \
| grep -P 'chr[0-9XY]+(?!_)'  > ${FILTERED_PEAK}






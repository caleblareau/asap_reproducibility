#!/bin/bash
#$ -S /bin/bash
#$ -cwd

SAMPLE_TMP=$1

mkdir -p output/genotype/pile_up_result output/genotype/final_res

echo SAMPLE : ${SAMPLE_TMP}

######### 
BAM_FILE=../../asap_large_data_files/CD4_stimulation_asapseq/cellranger/possorted_bam.bam;
BAM_FILE_2=output/genotype/chr1_22.bam;
VCF_FILE=output/genotype/C10.vcf;
BAM_FILE_3=output/genotype/filtered_for_popscle.bam;
BAM_FILE_3_SORT=output/genotype/filtered_for_popscle_sort.bam;
CellBarcode_LIST=output/genotype/cell_barcode_list.txt;
NUM_GENOTYPE=3

######### limit only autosome
samtools view -@ 30 -b ${BAM_FILE} chr{1..22} > ${BAM_FILE_2};
samtools index ${BAM_FILE_2};

######### freebayse
###### ref) https://github.com/ekg/freebayes.git
freebayes \
 -f ../../asap_large_data_files/CD4_stimulation_asapseq/reference/hg38.fa \
 -C 10 \
 ${BAM_FILE_2} > ${VCF_FILE}

bgzip ${VCF_FILE}
tabix -p vcf ${VCF_FILE}.gz

echo freebayes finished

######### making of cell-barcode list

cut -f 1,2 output/Signac/after_filter_Signac/HTO_res_filtered.txt \
| tail -n +2 > ${CellBarcode_LIST}

######### remake of bam file ## If it failed as background jobs, please do it as forground jobs.
######### ref) https://github.com/aertslab/popscle_helper_tools
filter_bam_file_for_popscle_dsc_pileup.sh \
 ${BAM_FILE_2} \
 ${CellBarcode_LIST} \
 ${VCF_FILE}.gz \
 ${BAM_FILE_3}

samtools index ${BAM_FILE_3};

######### popslce pileup
######### ref) https://github.com/statgen/popscle.git
popscle dsc-pileup \
 --sam ${BAM_FILE_3} \
 --vcf ${VCF_FILE}.gz \
 --out output/genotype/pile_up_result/result.pileup \
 --group-list ${CellBarcode_LIST}

echo Pile_up finished

######### freeemuxlet
~/anaconda3/envs/200624/bin/popscle freemuxlet \
 --plp output/genotype/pile_up_result/result.pileup \
 --out output/genotype/final_res/result \
 --nsample ${NUM_GENOTYPE}

echo freemuxlet finished 

zcat output/genotype/final_res/result.clust1.samples.gz \
| cut -f 2,13 > output/genotype/final_res/CellBarcode_genotype_result.txt;

rm ${BAM_FILE_2} ${BAM_FILE_2}.bai

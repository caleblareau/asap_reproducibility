#!/bin/bash
#$ -S /bin/bash

JOB_ID=$1;

WD=$(pwd);
SAMPLE_TMP=Perturb_CD4_stim;

echo SAMPLE : ${SAMPLE_TMP};
echo JOB_ID : ${JOB_ID};

mkdir -p output/GeneScore/bam/${JOB_ID};

######### dividing bam file into each CellBarcode
#### ref) https://github.com/timoast/sinto.git

cd ${WD}/output/GeneScore/bam/${JOB_ID} ;
sinto filterbarcodes  -b ${WD}/../../asap_large_data_files/CD4_stimulation_asapseq/cellranger/possorted_bam.bam \
                      -c ${WD}/output/GeneScore/job_list/cell_barcode_list_split${JOB_ID} \
                      -p 5 

cd ${WD}

######### bam -> sort -> BEDPE -> BED

for FILE_TMP in $(ls output/GeneScore/bam/${JOB_ID});do
 sgRNA_TMP=$(echo ${FILE_TMP} | cut -f 1,2 -d "_" -)
 CellBarcode_TMP=$(echo ${FILE_TMP} | cut -f 3 -d "_" - | cut -f 1 -d "-" -)

 BAM=output/GeneScore/bam/${sgRNA_TMP}/${FILE_TMP};
 BAM_SORT=output/GeneScore/tmp_bam/${sgRNA_TMP}/${FILE_TMP};
 BEDPE=output/GeneScore/bedpe/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}_v1.bedpe;

 TMP_BEDPE=output/GeneScore/tmp_bedpe/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}_v1_tmp.bedpe;
 TMP_BED=output/GeneScore/tmp_bed/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}_v1_tmp.bed;
 TMP_BED_2=output/GeneScore/tmp_bed/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}_v2_tmp.bed;
 BED=output/GeneScore/bed/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}.bed;
 BED_target=output/GeneScore/bed_target/${sgRNA_TMP}/${sgRNA_TMP}_${CellBarcode_TMP}_target.bed

 ##################
 BLACK_LIST=data/hg38-blacklist.v2.bed;
 TARGET_REGION_BED=data/GeneScore/gene_signature_position_list.bed
 ##################

 mkdir -p output/GeneScore/bam/${sgRNA_TMP} output/GeneScore/tmp_bam/${sgRNA_TMP};
 mkdir -p output/GeneScore/bedpe/${sgRNA_TMP} output/GeneScore/tmp_bedpe/${sgRNA_TMP};
 mkdir -p output/GeneScore/tmp_bed/${sgRNA_TMP} output/GeneScore/bed/${sgRNA_TMP};
 mkdir -p output/GeneScore/bed_target/${sgRNA_TMP};

 mv output/GeneScore/bam/${JOB_ID}/${FILE_TMP} ${BAM}

 samtools sort -n ${BAM} -o ${BAM_SORT};
 bedtools bamtobed -i ${BAM_SORT} -bedpe > ${BEDPE};

 awk 'OFS="\t"''{if($9=="+" && $10="-"){print $0}}' ${BEDPE} > ${TMP_BEDPE};
 awk 'OFS="\t"''{if($9=="-" && $10="+" && $1==$4){print $0}}' ${BEDPE} >> ${TMP_BEDPE};
 cut -f 1,2,6,7 ${TMP_BEDPE} \
 | awk 'OFS="\t"''{if($3<$2+2000 && $2<$3){print $0}}' -  > ${TMP_BED};

 for CHR_TMP in $(seq 1 22) X Y ; do
   grep -w chr${CHR_TMP} ${TMP_BED} >> ${TMP_BED_2}
 done

 sort -k 1,1 -k 2n,2 ${TMP_BED_2} \
 | bedtools intersect -v -a - -b ${BLACK_LIST} > ${BED}

 bedtools intersect \
    -a ${TARGET_REGION_BED} \
    -b ${BED} \
    -wa -c -sorted > ${BED_target}

 rm -r ${BAM} ${BAM_SORT} ${BEDPE} ${TMP_BEDPE} ${TMP_BED} ${TMP_BED_2}
done

rm -r output/GeneScore/bam/${JOB_ID}


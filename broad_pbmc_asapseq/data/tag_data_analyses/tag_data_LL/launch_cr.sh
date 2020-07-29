sample="1_NextGEM_PBMC_nodig"
bsub -q big-multi -n 12 /apps/lab/aryee/cellranger-atac-1.2.0/cellranger-atac count --fastqs NextGEM_PBMC_opt --id "${sample}_v12-mtMask" --sample "${sample}" --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0-mtMask --localcores 12 --force-cells 6000

sample="2_NextGEM_PBMC_plusdig"
bsub -q big-multi -n 12 /apps/lab/aryee/cellranger-atac-1.2.0/cellranger-atac count --fastqs NextGEM_PBMC_opt --id "${sample}_v12-mtMask" --sample "${sample}" --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0-mtMask --localcores 12 --force-cells 6000

sample="3_NextGEM_PBMC_nodigext"
bsub -q big-multi -n 12 /apps/lab/aryee/cellranger-atac-1.2.0/cellranger-atac count --fastqs NextGEM_PBMC_opt --id "${sample}_v12-mtMask" --sample "${sample}" --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0-mtMask --localcores 12 --force-cells 6000

sample="4_NextGEM_PBMC_plusdigext"
bsub -q big-multi -n 12 /apps/lab/aryee/cellranger-atac-1.2.0/cellranger-atac count --fastqs NextGEM_PBMC_opt --id "${sample}_v12-mtMask" --sample "${sample}" --reference /data/aryee/pub/genomes/cellranger/refdata-cellranger-atac-hg19-1.0.0-mtMask --localcores 12 --force-cells 6000




dir="1_NextGEM_PBMC_nodig_v12-mtMask"
bsub -q big-multi -n 12 mgatk bcall -i "${dir}/outs/possorted_bam.bam" -n "${dir}_mgatk" -o "${dir}_mgatk" -bt CB -b "${dir}/outs/filtered_peak_bc_matrix/barcodes.tsv" -c 12 -g rCRS

dir="2_NextGEM_PBMC_plusdig_v12-mtMask"
bsub -q big-multi -n 12 mgatk bcall -i "${dir}/outs/possorted_bam.bam" -n "${dir}_mgatk" -o "${dir}_mgatk" -bt CB -b "${dir}/outs/filtered_peak_bc_matrix/barcodes.tsv" -c 12 -g rCRS

dir="3_NextGEM_PBMC_nodigext_v12-mtMask"
bsub -q big-multi -n 12 mgatk bcall -i "${dir}/outs/possorted_bam.bam" -n "${dir}_mgatk" -o "${dir}_mgatk" -bt CB -b "${dir}/outs/filtered_peak_bc_matrix/barcodes.tsv" -c 12 -g rCRS

dir="4_NextGEM_PBMC_plusdigext_v12-mtMask"
bsub -q big-multi -n 12 mgatk bcall -i "${dir}/outs/possorted_bam.bam" -n "${dir}_mgatk" -o "${dir}_mgatk" -bt CB -b "${dir}/outs/filtered_peak_bc_matrix/barcodes.tsv" -c 12 -g rCRS





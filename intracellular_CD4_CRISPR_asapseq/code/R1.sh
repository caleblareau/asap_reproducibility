#!/bin/sh
#$ -S /bin/sh

code_R=$1
var1=$2

export LD_PRELOAD=/home/ytakeshima/.local/lib/libhdf5.so:/home/ytakeshima/.local/lib/libhdf5_hl.so
~/tool/R-3.6.3/bin/Rscript   --vanilla  --slave ${code_R} ${var1}


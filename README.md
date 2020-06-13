# ASAPseq paper reproducibility
Reproducing all analyses and figures for the ASAP-seq paper



## Reformatting fastqs for process with kite
We use kite (part of kallisto|bustools) to do all of the ADT/HTO tag quantification. 

Since the fastqs the come off the sequencer aren't currently supported by the kite
workflow, we have this accessory script to convert fastq files into a format similar 
to the 10x 5p scRNA-seq. This way, these can be more rapdily quantified. 

https://github.com/caleblareau/asap_to_kite



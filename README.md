# ASAPseq paper reproducibility
Reproducing all analyses and figures for the ASAP-seq paper


# Figure to directory mapping
**Note** this is for the updated version and _not_ the bioRxiv verison.
```
Figure 1 | species_mix_asapseq, pbmc_TBNK_comparisons_asapseq
Figure 2 | pbmc_TBNK_comparisons_asapseq
Figure 3 | bonemarow_asapseq
Figure 4 | pbmc_stimulation_asapseq_citeseq
Figure 5 | TBD
Figure 6 | CD4_CRISPR_asapseq
Figure 7 | intracellular_asapseq
```

## Reformatting fastqs for process with kite
We use kite (part of kallisto|bustools) to do all of the ADT/HTO tag quantification. 

Since the fastqs the come off the sequencer aren't currently supported by the kite
workflow, we have this accessory script to convert fastq files into a format similar 
to the 10x 5p scRNA-seq. This way, these can be more rapdily quantified. 

```
https://github.com/caleblareau/asap_to_kite
```



<br><br>
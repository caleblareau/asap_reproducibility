# Reproducing ASAP-seq and DOGMA-seq analyses

_Last updated: June 10, 2021_

This repository contains all code needed to reproduce all analyses:

**Scalable, multimodal profiling of chromatin accessibility, gene expression, and protein levels in single cells**
 
Eleni P. Mimitou+, Caleb A. Lareau+, Kelvin Y. Chen+ et al. Scalable, multimodal profiling of chromatin accessibility, gene expression, and protein levels in single cells. _Nature Biotechnology_ 2021. [Online here](https://www.nature.com/articles/s41587-021-00927-2)

# Setup

To best use this resource, we recommend pairing with large data files (that are not compatible with github as they exceed 100Mb). These files are available from the [Open Science Framework](https://osf.io/96va3/).

Once one downloads the zip archieve from OSF (~48 Gb), place the extracted folder in the same directory as this repository named `asap_large_data_files` (as shown below). This will enable running custom code to reproduce items in the output folders. 

```
.
├── asap_large_data_files
│   ├── CD4_stimulation_asapseq
│   ├── bonemarrow_data
│   ├── broad_experiment_pbmcs
│   ├── intraCD4_stimulation_asapseq
│   ├── multiome_pbmc_stim
│   ├── nygc_pbmc
│   ├── pbmc_stim_data
│   └── species_mix
└── asap_reproducibility
    ├── CD4_CRISPR_asapseq
    ├── README.md
    ├── bonemarow_asapseq
    ├── global_functions
    ├── intracellular_CD4_CRISPR_asapseq
    ├── intracellular_TBNK_asapseq
    ├── pbmc_TBNK_comparisons_asapseq
    ├── pbmc_stim_multiome
    ├── pbmc_stimulation_asapseq_citeseq
    └── species_mix_asapseq
```

The code assumes a relative file path with this organization. 


# Figure to directory mapping
**Note** this is for the version of the manuscript published in _Nature Biotechnolgy_ and _not_ the previous bioRxiv verison.
```
Figure 1 | species_mix_asapseq,pbmc_TBNK_comparisons_asapseq
Figure 2 | pbmc_TBNK_comparisons_asapseq
Figure 3 | bonemarow_asapseq
Figure 4 | pbmc_stimulation_asapseq_citeseq
Figure 5 | pbmc_stim_multiome
Figure 6 | CD4_CRISPR_asapseq
Figure 7 | intracellular_asapseq,intracellular_CD4_CRISPR_asapseq
```

## Reformatting fastqs for process with kite
We use kite (part of kallisto|bustools) to do all of the ADT/HTO tag quantification. 

Since the fastqs the come off the sequencer aren't currently supported by the kite
workflow, we have this accessory script to convert fastq files into a format similar 
to the 10x 5p scRNA-seq. This way, these can be more rapdily quantified. 

```
https://github.com/caleblareau/asap_to_kite
```

## Questions/ contact

Formal questions about this work should be addressed to the corresponding author: [Peter Smibert](mailto:smibertp@gmail.com). 

For questions or comments related to the specific code, please raise an issue or submit a pull request. 

<br><br>

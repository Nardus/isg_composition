# Prediction of interferon stimulated and repressed genes using sequence composition

[![DOI](https://zenodo.org/badge/378905292.svg)](https://zenodo.org/badge/latestdoi/378905292)


Code for the machine learning component of Shaw, Rihn, Mollentze, *et al.* (2021) "The ‘antiviral state’ has shaped the CpG composition of the vertebrate interferome". Raw data for this publication are available from DOIs [10.5525/gla.researchdata.1159](http://doi.org/10.5525/gla.researchdata.1159) and [10.5281/zenodo.5035606](https://doi.org/10.5281/zenodo.5035606), and will be downloaded automatically as needed by this analysis pipeline. 

The aim of this part of the analysis was to investigate to what extent interferon-stimulated genes, interferon-repressed genes, and remaining genes are compositionally distinct, and to identify the most important compositional features allowing them to be distinguished. 


## Usage
Analyses were run using R version 3.5.1. Required R libraries are managed by packrat. A full run requires approximately 15GB of disk space, and by default will run in parallel across 8 threads.

To re-run the entire pipeline, use 
```
Rscript packrat/init.R --args --bootstrap-packrat
make all
```

Run `make help` for further details.



## Folder structure
`[*]` Indicates folders which will be created while running the pipeline

```
└─isg_composition/
   ├─CalculatedData/...................... [*] Final pre-processed data used to train models
   │ 
   ├─Data/ ............................... Raw data (most files downloaded from the main data  
   │   │                                   repositories as needed)
   │   ├─ISG_identities.csv                Interferome data, from Shaw et al. 2017 (downloaded
   │   │                                   from https://isg.data.cvr.ac.uk)
   │   ├─consistent_irgs.csv               IRGs consistently found across all experiments
   │   ├─MouseHoldout/                     [*] Holdout data describing the mouse interferome
   │   │                                       (from Dölken et al., 2008)
   │   ├─Other_Experiments/                [*] Data from the knockdown- and other experiments
   │   │                                       described in the manuscript
   │   ├─Part1_Dinucs/                     [*] Expression data and dinucleotide biases (see 
   │   │                                       "Data" below)
   │   ├─Part2_Codons/                     [*] Codon and codon-pair biases, etc
   │   └─SelectedGenes/                    [*] Genes analyzed in the manuscript
   │
   ├─Output/ ............................. [*] Generated output files
   │
   ├─packrat/ ............................ Record of R libraries used
   │
   ├─Scripts/ ............................ Main pipeline scripts
   │
   └─Makefile ............................ Record of workflow and dependencies between files
```

### Data
Data will be downloaded as needed while running the pipeline, but can also be obtained using `make download_data` (but note that not all data from the manuscript will be extracted, and most files in the "Other_Experiments" folder will be renamed - refer to the original repositories linked above for the canonical versions of these files). The core genome composition data retreived are split across several files:

- In `Data/Part1_Dinucs/`:
	- This folder contains the expression data and dinucleotide composition summaries for the **longest transcript** of each gene
  - Files named `[species]_cds_new_dat_fpkm_dups.txt` contains data for coding sequences only
  - Files named `[species]_cdna_new_dat_fpkm_dups.txt` contains the same calculations across the entire coding sequence
  - Most files match cds/cdna datasets provided by [https://www.ensembl.org/info/data/ftp/index.html]
    - *Where [species] is replaced with 'biomart':* For humans only, contains **all** human genes, downloaded manually from Ensembl Biomart
- In `Data/Part2_Codons/`:
  - Contains codon and codon-pair biases for the same transcripts as above
  - File names have the format `[species]_cds_new_cpb_dat_dups.txt`
  - Note that there is no corresponding `cdna` version here, since these features are valid for coding regions only

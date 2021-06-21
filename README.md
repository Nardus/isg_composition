# Prediction of interferon stimulated and repressed genes using sequence composition

Code for the machine learning component of Shaw, Rihn, Mollentze, *et al.* (2021) "The ‘antiviral state’ has shaped the CpG composition of the vertebrate interferome". Raw data for this publication are available from ___, and will be downloaded automatically by this analysis pipeline. 

The aim of this part of the analysis was to investigate to what extent interferon-stimulated genes, interferon-repressed genes, and remaining genes are compositionally distinct, and to identify the most important compositional features allowing them to be distinguished. 


## Usage
Analyses were run using R version 3.5.1. Required R libraries are managed by packrat.

To re-run the entire pipeline, use 
```
Rscript -e "packrat::restore()"
make all
```

Run `make help` for further details.


## Folder structure
`[*]` Indicates folders which will be created while running the pipeline

```
└─isg_composition/
   ├─CalculatedData/................................ [*] Final pre-processed data used to train models
   │ 
   ├─Data/ ......................................... Raw data (downloaded as needed)
   │   │
   │   ├─A549Holdout/
   │   ├─MouseHoldout/
   │   ├─Other_Experiments/
   │   ├─Part1_Dinucs/                               Expression data and dinucleotide biases 
   │   │                                             (see "Composition data" below)
   │   ├─Part2_Codons/                               Codon and codon-pair biases, etc
   │   ├─SelectedGenes/                              Genes selected analyzed in the manuscript
   │   └─TranscriptMetadata/
   │
   ├─Output/ ....................................... [*] Generated output files
   │   │
   │   ├─Plots/                                      Result and diagnostic plots, including panels
   │   │                                             for figure 1
   │   └─predictions/                                Model predictions displayed in other figures
   │
   ├─packrat/ ......................................
   │
   ├─Scripts/ ...................................... Main pipeline scripts
   │
   └─Makefile ...................................... Record of workflow and dependencies between files
```

### Composition data
After running `make download_data`, sevaral genome composition files will be available

- In `Data/Part1_Dinucs/`:
		- This folder contains the expression data and several dinucleotide summaries for the **longest transcript** of each gene
		    - For all 16 dinucleotides:
		        - Observed/Expected ratio - labelled e.g. ApA etc;
		        - Raw dinucleotide percentages - labelled e.g. AA%
    - Files named `[host]_cds_new_dat_fpkm_dups.txt` contains data for coding sequences only
    - Files named `[host]_cdna_new_dat_fpkm_dups.txt` contains the same calculations across the entire coding sequence
    - Most files match cds/cdna datasets provided by [https://www.ensembl.org/info/data/ftp/index.html]
        - *Where [host] is replaced with 'biomart':* For humans only, contains **all** human genes, downloaded manually from Ensembl Biomart
- In `Data/Part2_Codons/`:
    - Contain summaries of codon and codon-pair biases, etc.
    - File names have the format `[host]_cds_new_cpb_dat_dups.txt` as in Part1
    - Note that there is no corresponding `cdna` version here, since these features are valid for coding regions only

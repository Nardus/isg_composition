# Prediction of interferon stimulated and repressed genes using sequence composition

Code for the machine learning component of Shaw, Rihn, Mollentze, *et al.* (2021) "The ‘antiviral state’ has shaped the CpG composition of the vertebrate interferome". 

The aim of this part of the analysis was to investigate to what extent interferon-stimulated genes, interferon-repressed genes, and remaining genes are compositionally distinct, and to identify the most important compositional features allowing them to be distinguished. The raw data needed to run this pipeline are available separately (see manuscript for details). 


## Usage
Analyses were run using R version 3.5.1. Required R libraries are managed by packrat.

To re-run the entire pipeline, use 
```
Rscript -e "packrat::restore()"
make all
```

Run `make help` for further details.


## Folder structure
- `[*]` Indicates folders which will be created while running the pipeline
- `[**]` The data required to run this pipeline are available separately (see manuscript for details)

```
└─isg_composition/
   ├─CalculatedData/...................... [*] Final pre-processed data used to train models
   │ 
   ├─Data/ ............................... [**] Place downloaded data here  
   │
   ├─Output/ ............................. [*] Generated output files
   │
   ├─packrat/ ............................ Record of R libraries used
   │
   ├─Scripts/ ............................ Main pipeline scripts
   │
   └─Makefile ............................ Record of workflow and dependencies between files
```

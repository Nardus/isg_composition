# -------------------------------------------------------------------------------------------------
#' Main pipeline for ISG prediction from their genome composition
#   - Data is assumed to exist already
# -------------------------------------------------------------------------------------------------
#'
#' This pipeline can run in four modes:
#'  1. Training on the top 50 genes from each class
#'     1.1. Across all hosts simultaneously (used in the manuscript)
#'     1.2. Genes from each host individually
#'  2. Training on the ALL genes from each class
#'     2.1. Across all hosts simultaneously
#'     2.2. Genes from each host individually
#'
#' In practice, results are largely similar, but using genes from each host  
#' individually (1.2 & 2.2) has smaller sample sizes during training. In the 
#' manuscript, option 1.1 was used to focus on the most extreme differences, 
#' which is diluted when including all genes (many of which are only marginally
#' overexpressed/repressed). This is described in more detail in the manuscript.
#'
#' Usage:
#'   "make all" performs the analyses described in the manuscript.
#'   Use "make additional", to run the alternative modes (currently implemented
#'   for specific combinations of genes/hosts only; edit the preamble of the
#'   Makefile if other combinations are needed).
#'

.NOTPARALLEL:  # Scripts below already parallelised, so force MAKE to run sequentially


# ---- Expected output ----------------------------------------------------------------------------
# (Change these to trigger additional runs)

# Hosts, rougly in order of increasing distance to humans (useful for plot order):
HOSTS = human rat megabat microbat cow pig sheep horse dog chicken

PREPARED_DATA_HOSTS = $(patsubst %, CalculatedData/InputData_%.rds, $(HOSTS))

## TRAINING - ALL HOSTS COMBINED (top 50 only):
TRAIN_HOSTS_COMBINED = Output/AllHostsCombined/ISGvsIRG-Top50/Fit_ISGvsIRG-Top50_AllHosts.rds \
						Output/AllHostsCombined/ISGvsRandom-Top50/Fit_ISGvsRandom-Top50_AllHosts.rds \
						Output/AllHostsCombined/IRGvsRandom-Top50/Fit_IRGvsRandom-Top50_AllHosts.rds \
						Output/AllHostsCombined/3Class-Top50/Fit_3Class-Top50_AllHosts.rds

## ADDITIONAL TRAINING - BY HOST:
# All genes: doing this for human only
TRAIN_ISGvsIRG_AllGenes = Output/ISGvsIRG-AllGenes/Fit_ISGvsIRG-AllGenes_human.rds
TRAIN_ISGvsRandom_AllGenes = Output/ISGvsRandom-AllGenes/Fit_ISGvsRandom-AllGenes_human.rds
TRAIN_IRGvsRandom_AllGenes = Output/IRGvsRandom-AllGenes/Fit_IRGvsRandom-AllGenes_human.rds

# Top 50:
TRAIN_ISGvsIRG_Top50 = $(patsubst %, Output/ISGvsIRG-Top50/Fit_ISGvsIRG-Top50_%.rds, $(HOSTS))
TRAIN_ISGvsRandom_Top50 = Output/ISGvsRandom-Top50/Fit_ISGvsRandom-Top50_human.rds  # Human-only for now
TRAIN_IRGvsRandom_Top50 = Output/IRGvsRandom-Top50/Fit_IRGvsRandom-Top50_human.rds


ADDITIONAL_RUNS = $(TRAIN_ISGvsIRG_AllGenes) $(TRAIN_ISGvsRandom_AllGenes) $(TRAIN_IRGvsRandom_AllGenes) \
			 	  $(TRAIN_ISGvsIRG_Top50) $(TRAIN_ISGvsRandom_Top50) $(TRAIN_IRGvsRandom_Top50)


# ---- Input files --------------------------------------------------------------------------------
IDENTITY_DATA = $(patsubst %, Data/SelectedGenes/%_Top50_Formatted.csv, $(HOSTS))

FEATURE_DATA = $(patsubst %, Data/Part1_Dinucs/%_cds_new_dat_fpkm_dups.txt, $(HOSTS)) \
			   $(patsubst %, Data/Part1_Dinucs/%_cdna_new_dat_fpkm_dups.txt, $(HOSTS)) \
			   $(patsubst %, Data/Part2_Codons/%_cds_new_cpb_dat_dups.txt, $(HOSTS))

MOUSE_HOLDOUT_DATA = Data/MouseHoldout/mouse_holdout_genes.csv \
					 Data/MouseHoldout/mouse_cds.txt \
					 Data/MouseHoldout/mouse_cds_cg.txt \
					 Data/MouseHoldout/mouse_cdna_ncrna.txt

EXPERIMENT_DATA = Data/Other_Experiments/consistent_irgs.csv \
				  Data/Other_Experiments/Random100_nonDE_CPM1_A549_typeI_4h.csv \
				  Data/Other_Experiments/EdgeRCas9ClonesIFNome_withFDR_forPrism.xlsx \
				  Data/Other_Experiments/EdgeR_KOClonesIFNome_withFDR_forPrisim.xlsx \
				  Data/Other_Experiments/Fig_5I_Top50-Up-Down-Random.csv \
				  Data/ISG_identities.csv \
				  Data/SelectedGenes/human_Top50_Formatted.csv


# ---- Commands -----------------------------------------------------------------------------------
## SHORTCUTS:
# Main analysis
.PHONY: all, download_data, prepare, train, predict, plot
all: train predict plot  # Not including "download_data" and "prepare" here: files should extracted/processed as needed for training instead

download_data: $(IDENTITY_DATA) $(FEATURE_DATA) $(MOUSE_HOLDOUT_DATA) $(EXPERIMENT_DATA)

prepare: $(PREPARED_DATA_HOSTS) CalculatedData/InputData_AllHosts.rds

train: $(TRAIN_HOSTS_COMBINED)

predict: Data/MouseHoldout/mouse_holdout_predictions.csv

plot: Output/AllHostsCombined/Plots/OverallAccuracy_Combined.pdf \
	  Output/AllHostsCombined/Plots/cpg+probability_ISGvsIRG_model.pdf \
	  Output/AllHostsCombined/Plots/cpg+probability_IRGvsRandom_model.pdf
	@echo "Done"


# Additional combinations
.PHONY: additional, train_additional, plot_additional
additional: train_additional, plot_additional

train_additional: $(ADDITIONAL_RUNS)

plot_additional: Output/Plots/OverallAccuracy_ISGvsIRG.pdf
	@echo "Done"


## ACTUAL COMMANDS:
# Note: the '$*' below expands to the part of the filename that was matched
#       by '%' in the rule, in this case the host name

# Download and extract data
Data/main_data.zip:
	curl -L -o $@ --create-dirs http://researchdata.gla.ac.uk/1159/1/Data_1159.zip

Data/other_data.zip:
	curl -L -o $@ --create-dirs https://zenodo.org/record/5035607/files/rjorton/DinucISG-v1.0.0.zip

# - Identity data (needs some renaming)
HOST_CAPITALIZED = $(shell echo $* | awk '{ print toupper(substr($$0,1,1)) substr($$0,2) }')

Data/SelectedGenes/%_Top50_Formatted.csv: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 2/Fig 2 data/${HOST_CAPITALIZED}_Top50_Formatted.csv" > $@

Data/SelectedGenes/megabat_Top50_Formatted.csv: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 2/Fig 2 data/Fruit bat_Top50_Formatted.csv" > $@


# - Part 1 feature files
Data/Part1_Dinucs/%_cdna_new_dat_fpkm_dups.txt: Data/other_data.zip
	mkdir -p $(@D)
	unzip -p $< rjorton-DinucISG-d67e82b/DinucData/$*_cdna_new_dat_fpkm_dups.txt > $@

Data/Part1_Dinucs/%_cds_new_dat_fpkm_dups.txt: Data/other_data.zip
	mkdir -p $(@D)
	unzip -p $< rjorton-DinucISG-d67e82b/DinucData/$*_cds_new_dat_fpkm_dups.txt > $@


# - Part 2 feature files
Data/Part2_Codons/%_cds_new_cpb_dat_dups.txt: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Sequence Features/Enlighten_ISG_Zap/$*_cds_new_cpb_dat_dups.txt" > $@


# - Mouse holdout data
Data/MouseHoldout/mouse_holdout_genes.csv: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 1/Fig 1H/Murine genes.csv" > $@

Data/MouseHoldout/mouse_cds.txt: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Sequence Features/Enlighten_ISG_Zap/mouse_cds.txt" > $@

Data/MouseHoldout/mouse_cds_cg.txt: Data/other_data.zip
	mkdir -p $(@D)
	unzip -p $< rjorton-DinucISG-d67e82b/DinucData/mouse_cds_cg.txt > $@

Data/MouseHoldout/mouse_cdna_ncrna.txt: Data/other_data.zip
	mkdir -p $(@D)
	unzip -p $< rjorton-DinucISG-d67e82b/DinucData/mouse_cdna_ncrna.txt > $@


# - Experiment data
Data/Other_Experiments/consistent_irgs.csv: Data/consistent_irgs.csv
	mkdir -p $(@D)
	cp $< $@

Data/Other_Experiments/Random100_nonDE_CPM1_A549_typeI_4h.csv: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 2/Fig 2 data/Random100_nonDE_CPM1_A549_typeI_4h.csv" > $@

Data/Other_Experiments/EdgeRCas9ClonesIFNome_withFDR_forPrism.xlsx: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 4/Fig 4DE/Fig 4DE No Guide Interferomes.xlsx" > $@

Data/Other_Experiments/EdgeR_KOClonesIFNome_withFDR_forPrisim.xlsx: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 4/Fig 4DE/Fig 4DE ZAP KO interferome.xlsx" > $@

Data/Other_Experiments/Fig_5I_Top50-Up-Down-Random.csv: Data/main_data.zip
	mkdir -p $(@D)
	unzip -p $< "ENLIGHTEN - Shaw, Rihn, Mollentze et al., 2021/Fig 5/Fig 5I/Fig 5I Top50-Up-Down-Random.csv" > $@



# Prepare data
CalculatedData/InputData_AllHosts.rds: $(FEATURE_DATA) Data/ISG_identities.csv
	mkdir -p $(@D)
	Rscript --default-packages=methods,utils,stats ./Scripts/Prepare.R --host $(HOSTS)

CalculatedData/InputData_%.rds: Data/Part1_Dinucs/%_cds_new_dat_fpkm_dups.txt \
								Data/Part1_Dinucs/%_cdna_new_dat_fpkm_dups.txt \
								Data/Part2_Codons/%_cds_new_cpb_dat_dups.txt \
								Data/ISG_identities.csv
	mkdir -p $(@D)
	Rscript --default-packages=methods,utils,stats ./Scripts/Prepare.R --host $*


# Training:
#  - All genes
Output/3Class-AllGenes/Fit_3Class-AllGenes_%.rds: CalculatedData/InputData_%.rds \
												  Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --Predict3Class --host $*

Output/ISGvsIRG-AllGenes/Fit_ISGvsIRG-AllGenes_%.rds: CalculatedData/InputData_%.rds \
													  Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsIRG --host $*


Output/ISGvsRandom-AllGenes/Fit_ISGvsRandom-AllGenes_%.rds: CalculatedData/InputData_%.rds \
															Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsRandom --host $*

Output/IRGvsRandom-AllGenes/Fit_IRGvsRandom-AllGenes_%.rds: CalculatedData/InputData_%.rds \
															Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --IRGvsRandom --host $*


#  - Top 50
Output/ISGvsIRG-Top50/Fit_ISGvsIRG-Top50_%.rds: CalculatedData/InputData_%.rds \
												Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsIRG --Top50 --host $*

Output/ISGvsRandom-Top50/Fit_ISGvsRandom-Top50_%.rds: CalculatedData/InputData_%.rds \
													  Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsRandom --Top50 --host $*

Output/IRGvsRandom-Top50/Fit_IRGvsRandom-Top50_%.rds: CalculatedData/InputData_%.rds \
													  Data/SelectedGenes/%_Top50_Formatted.csv
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --IRGvsRandom --Top50 --host $*


# - All hosts combined:
Output/AllHostsCombined/ISGvsIRG-Top50/Fit_ISGvsIRG-Top50_AllHosts.rds: CalculatedData/InputData_AllHosts.rds \
																		$(IDENTITY_DATA)
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsIRG --Top50 --host AllHosts
	mkdir -p Output/AllHostsCombined/ISGvsIRG-Top50
	mv Output/ISGvsIRG-Top50/Fit_ISGvsIRG-Top50_AllHosts.rds $(@D)/
	mv Output/ISGvsIRG-Top50/Data_ISGvsIRG-Top50_AllHosts.rds $(@D)/

Output/AllHostsCombined/ISGvsRandom-Top50/Fit_ISGvsRandom-Top50_AllHosts.rds: CalculatedData/InputData_AllHosts.rds \
																			  $(IDENTITY_DATA)
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --ISGvsRandom --Top50 --host AllHosts
	mkdir -p Output/AllHostsCombined/ISGvsRandom-Top50
	mv Output/ISGvsRandom-Top50/Fit_ISGvsRandom-Top50_AllHosts.rds $(@D)/
	mv Output/ISGvsRandom-Top50/Data_ISGvsRandom-Top50_AllHosts.rds $(@D)/

Output/AllHostsCombined/IRGvsRandom-Top50/Fit_IRGvsRandom-Top50_AllHosts.rds: CalculatedData/InputData_AllHosts.rds \
																			  $(IDENTITY_DATA)
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --IRGvsRandom --Top50 --host AllHosts
	mkdir -p Output/AllHostsCombined/IRGvsRandom-Top50
	mv Output/IRGvsRandom-Top50/Fit_IRGvsRandom-Top50_AllHosts.rds $(@D)/
	mv Output/IRGvsRandom-Top50/Data_IRGvsRandom-Top50_AllHosts.rds $(@D)/

Output/AllHostsCombined/3Class-Top50/Fit_3Class-Top50_AllHosts.rds: CalculatedData/InputData_AllHosts.rds \
																	$(IDENTITY_DATA)
	Rscript --default-packages=methods,utils,stats Scripts/Train.R --Predict3Class --Top50 --host AllHosts
	mkdir -p Output/AllHostsCombined/3Class-Top50
	mv Output/3Class-Top50/Fit_3Class-Top50_AllHosts.rds $(@D)/
	mv Output/3Class-Top50/Data_3Class-Top50_AllHosts.rds $(@D)/


# Predict:
Data/MouseHoldout/mouse_holdout_predictions.csv: $(MOUSE_HOLDOUT_DATA)
	Rscript Scripts/PredictMouseHoldout.R


# Plots:
# - Main (all hosts combined)
#   - All diagnostic plots
Output/AllHostsCombined/Plots/OverallAccuracy_Combined.pdf: $(TRAIN_HOSTS_COMBINED)
	Rscript --default-packages=methods,utils,stats,grDevices ./Scripts/MakePlots_MultiHost.R

#   - Predictions from A549 experiments (figure 4/5)
Output/AllHostsCombined/Plots/cpg+probability_ISGvsIRG_model.pdf: $(TRAIN_HOSTS_COMBINED) \
																  $(EXPERIMENT_DATA) \
																  CalculatedData/InputData_biomart.rds
	Rscript Scripts/MakePlots_KnockdownExperiments.R ISGvsIRG ISG

Output/AllHostsCombined/Plots/cpg+probability_IRGvsRandom_model.pdf: $(TRAIN_HOSTS_COMBINED) \
																	 $(EXPERIMENT_DATA) \
																	 CalculatedData/InputData_biomart.rds
	Rscript Scripts/MakePlots_KnockdownExperiments.R IRGvsRandom IRG


# - Additional runs (separate hosts, etc.)
Output/Plots/OverallAccuracy_ISGvsIRG.pdf: $(ADDITIONAL_RUNS)
	Rscript --default-packages=methods,utils,stats,grDevices ./Scripts/MakePlots.R $(HOSTS);



# ---- Auto document this file --------------------------------------------------------------------
# Comments above that start with #' become help strings

.PHONY: help
help: Makefile
	@grep "^#'" $< | cut -c4-

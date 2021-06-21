#!Rscript

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##   ISGpredict - Extra
##     Predict ISGs for a holdout species / independent dataset
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(rprojroot)
library(tidyverse)
library(xgboost)
library(caret)


# Data column formats:
NUC_CONTENT <- '^(GC|AT)cont$'
DINUC_BIAS <- '^[ATGCU]p[ATGCU]$'
DINUC_PERCENT <- '^[ATGCU]{2}%$'

BRIDGE_DINUC <- '^br[ATGCU]p[ATGCU]$'
NON_BRIDGE_DINUC <- '^NonBr[ATGCU]p[ATGCU]$'
#DINUC_RAW <- '^Raw[ATGCU]{2}$'      # These seem to be the same as DINUC_PERCENT in part1
AMINO_BIAS <- '^[[:alpha:]]-Bias$'
CODON_BIAS <- '^[ATGCU]{3}-Bias$'

CODON_PAIR <- '^[ATCGU]{3}\\([[:alpha:]]\\)-[ATCGU]{3}\\([[:alpha:]]\\)$'  # Not currently used


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Paths relative to project root:
rootDir <- find_rstudio_root_file()
UtilsFile <- file.path(rootDir, 'Scripts', 'Utils.R')
holdoutDir <- file.path(rootDir, 'Data', 'MouseHoldout')
metaDataPath <- file.path(holdoutDir, 'mouse_holdout_genes.csv')

source(UtilsFile)

# Read metadata:
metaData <- read_csv(metaDataPath, na = c('NA', '-'))

# Read feature data:
#  - mouse_cds.txt: contains all features, calculated for coding sequences
#  - mouse_cdna_ncrna.txt: contains the same features (where applicable), calculated across entire cdna. 
# 	 This also means ncrna genes can be included here
codingFeatures <- read_data('mouse', suffix = '_cds.txt', directory = holdoutDir) %>% 
	distinct()

codingExtra <- read_data('mouse', suffix = '_cds_cg.txt', directory = holdoutDir) %>% 
	distinct()

geneFeatures <- read_data('mouse', suffix = '_cdna_ncrna.txt', directory = holdoutDir) %>% 
	distinct()


# Read trained model:
modelPath <- file.path(rootDir, 'Output', 'AllHostsCombined', 'ISGvsIRG-Top50', 'Fit_ISGvsIRG-Top50_AllHosts.rds')
trainedModel <- readRDS(modelPath)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Find top 50 ISGs/IRGs ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Some genes have multiple entries - we want the top 50 *unique* genes
topISG <- metaData %>% 
	group_by(`Ensembl ID`) %>% 
	summarise(MaxChange = max(`Fold Change`)) %>% 
	ungroup() %>% 
	top_n(MaxChange, n = 50) %>% 
	mutate(Label = 'Stimulated') %>% 
	select(-MaxChange)

topIRG <- metaData %>% 
	group_by(`Ensembl ID`) %>% 
	summarise(MinChange = max(`Fold Change`)) %>% 
	ungroup() %>% 
	top_n(MinChange, n = -50) %>% 
	mutate(Label = 'Repressed') %>% 
	select(-MinChange)

stopifnot(!any(topISG$`Ensembl ID` %in% topIRG$`Ensembl ID`))


metaData <- bind_rows(topISG, topIRG)
stopifnot(length(unique(metaData$`Ensembl ID`)) == nrow(metaData))


write_csv(metaData, path = file.path(holdoutDir, 'mouse_holdout_genes_top50unique.csv'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Rename columns to id data source -----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Columns to keep / rename
nameColumns <- c('^Gene$') %>% 
	paste(collapse = '|')

codingColumns <- c(nameColumns, DINUC_BIAS, BRIDGE_DINUC, NON_BRIDGE_DINUC, AMINO_BIAS, CODON_BIAS) %>% 
	paste(collapse = '|')

codingExtraColumns <- c(nameColumns, NUC_CONTENT, DINUC_PERCENT) %>% 
	paste(collapse = '|')

geneColumns <- c(nameColumns, NUC_CONTENT, DINUC_BIAS, DINUC_PERCENT) %>% 
	paste(collapse = '|')


# Rename:
codingFeatures <- codingFeatures %>% 
	mutate(SeqName = str_remove(SeqName, '^>')) %>% 
	separate(SeqName, into = c('Gene', 'Transcript'), sep = '\\|') %>% 
	select(matches(codingColumns, ignore.case = FALSE)) %>% 
	rename(ApT = ApU,
				 CpT = CpU,
				 GpT = GpU,
				 UpT = UpU) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_Coding'))

codingExtra <- codingExtra %>% 
	select(matches(codingExtraColumns, ignore.case = FALSE)) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_Coding'))
	
geneFeatures <- geneFeatures %>% 
	select(matches(geneColumns, ignore.case = FALSE)) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_EntireGene'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Merge data ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
joinedData <- codingFeatures %>% 
	full_join(codingExtra, by = 'Gene') %>% 
	full_join(geneFeatures, by = 'Gene')

# Check for duplication:
stopifnot(length(unique(joinedData$Gene)) == nrow(joinedData))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Predict ------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
metaData <- metaData %>% 
	distinct(`Ensembl ID`, Label)


preds <- data.frame(Gene = joinedData$Gene,
										PredictedLabel = predict(trainedModel, newdata = joinedData, na.action = 'na.pass')) %>% 
	left_join(metaData, by = c('Gene' = 'Ensembl ID'))


# Check accuracy:
preds %>% 
	group_by(Label) %>% 
	summarise(Accuracy = sum(PredictedLabel == Label),
						Proportion_accurate = Accuracy / n())

# Save
preds %>% 
	select(Gene, ActualLabel = Label, PredictedLabel) %>% 
	write_excel_csv(file.path(holdoutDir, 'mouse_holdout_predictions.csv'))

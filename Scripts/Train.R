#!Rscript

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##   ISGpredict - Part 2
##     Tune and train models
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(rprojroot)
library(tidyverse)
library(xgboost)
library(caret)
library(doParallel)
library(argparse)


# Read user input:
parser <- ArgumentParser()
parser$add_argument('--host', 
										help = 'Host to use in training. Assumes Prepare.R has already been called for this host. The words AllHosts and human trigger special handling when identifying the top 50 genes.')

parser$add_argument('--Predict3Class', action = "store_true", default = FALSE,
										help = 'Predict ISG vs IRG vs Random. Can be combined with other flags below to fit multiple models')
parser$add_argument('--ISGvsIRG', action = "store_true", default = FALSE,
										help = 'Predict ISG vs IRG.')
parser$add_argument('--ISGvsRandom', action = "store_true", default = FALSE,
										help = 'Predict ISG vs random (non DE) genes')
parser$add_argument('--IRGvsRandom', action = "store_true", default = FALSE,
										help = 'Predict IRG vs random (non DE) genes')

parser$add_argument('--Top50', action = "store_true", default = FALSE,
										help = 'Use top 50 genes from each set only (default: use all genes)')


INPUT_ARGS <- parser$parse_args()

if (!(INPUT_ARGS$Predict3Class | INPUT_ARGS$ISGvsIRG | INPUT_ARGS$ISGvsRandom | INPUT_ARGS$IRGvsRandom))
	stop('Nothing to predict. Select at least one flag (use Train.R --help to see valid options)')


# OTHER CONSTANTS:
# Repeated cross-validation settings:
CV_FOLDS <- 5    # Number of cross-validation folds
CV_TIMES <- 50     # Number of repeated CV's to perform

# Parallel settings:
N_CORES <- 8


TOP_DIR <- find_rstudio_root_file()

OUT_DIR <- file.path(TOP_DIR, 'Output')
dir.create(OUT_DIR, recursive = T)

cl <- makeCluster(N_CORES)
registerDoParallel(cl)


load_top_genes <- function(host) {
	# Load top genes for a single host
	topfile <- paste0(host, '_Top50_Formatted.csv') %>% 
		file.path('Data', 'SelectedGenes', .)
	
	read_csv(topfile, col_types = cols(.default = 'c')) %>% 
		select(Gene = `Ensembl gene ID`, ExpressionCategory = Class) %>% 
		mutate(HOST = host,
					 ExpressionCategory = if_else(.data$ExpressionCategory == 'Random', 'NonDE', 
					 														  if_else(.data$ExpressionCategory == 'ISG', 'Stimulated', 
					 														  				'Repressed')))
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
dataFile <- paste0('InputData_', INPUT_ARGS$host, '.rds')
dataPath <- file.path(TOP_DIR, 'CalculatedData', dataFile)

inputData <- readRDS(dataPath)

if (length(unique(inputData$HOST)) != 1) {
	# Multiple hosts is valid, but need to make sure it was what we intended to do,
	# in part because the model will not be given this information (HOST column 
	# removed below)
	unique(inputData$HOST) %>% 
		paste(collapse = ', ') %>% 
		paste('Input contains multiple hosts:', .) %>% 
		warning()
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Setup --------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
response <- 'Label'
nonFeatures <- c('HOST', 'Gene', 'Log2FC', 'FDR')

trainSettings <- trainControl(method = 'repeatedcv',
															number = CV_FOLDS,
															repeats = CV_TIMES,
															search = 'random',
															returnResamp = 'all',
															sampling = 'down',
															savePredictions = 'final',
															classProbs = TRUE,
															summaryFunction = mnLogLoss)

# Primary data consists of all genes observed to be expressed:
if (! INPUT_ARGS$Top50) {
	OutNameBase <- 'AllGenes'
	
	FinalData <- inputData %>% 
		filter(!is.na(Log2FC) & !is.na(ExpressionCategory))
	
} else {
	# Reduce this to top 50:
	OutNameBase <- 'Top50'
	
	topGenes <- unique(inputData$HOST) %>% 
		lapply(FUN = load_top_genes) %>% 
		bind_rows()
	
	FinalData <- inputData %>% 
		select(-.data$ExpressionCategory) %>% 
		inner_join(topGenes, by = c('HOST', 'Gene'))
	
	stopifnot(nrow(FinalData) == 50*3*length(unique(inputData$HOST)))  # Expect three categories of 50 genes for each host
}


# Diagnostic message:
message(sprintf('Final input dataset contains %d observations from %d species', 
								nrow(FinalData), 
								length(unique(FinalData$HOST))))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Training -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
save_object <- function(obj, runid, type = 'Data') {
	# Save data/fit using the same basic pattern for all filenames
	outName <- paste(type, runid, INPUT_ARGS$host, sep = '_')  # e.g. Data_ISGvsIRG_human
	outDir <- file.path(OUT_DIR, runid)
	dir.create(outDir, recursive = T, showWarnings = F)
	saveRDS(obj, file = file.path(outDir, paste0(outName, '.rds')))
}

train_model <- function(data, runid, removeCols = nonFeatures) {
	# Train model and save fit to disk
	message(sprintf('Training model on %d observations from %d species', 
									nrow(data), 
									length(unique(data$HOST))))
	
	data <- data %>% 
		select(-one_of(removeCols))
	
	trainedModel <- train(ExpressionCategory ~ .,
												data = data,
												method = 'xgbDART', 
												trControl = trainSettings,
												metric = 'logLoss',
												tuneLength = 100,
												na.action = 'na.pass',
												nthread = 1)
	
	save_object(trainedModel, runid = runid, type = 'Fit')
}


## ---- 3-way prediction ---------------------------------------------------------------------------
if (INPUT_ARGS$Predict3Class) {
	# This uses all data
	# - Try to distinguish all three possible class (up/down/no change)
	classData <- FinalData
	classData$ExpressionCategory <- factor(classData$ExpressionCategory, levels = c('NonDE', 'Stimulated', 'Repressed')) # Make NonDE the baseline category
	
	runid <- paste(c('3Class', OutNameBase), collapse = '-')
	save_object(classData, runid = runid, type = 'Data')
	
	# Training:
	train_model(classData, runid = runid)
}



## ---- ISG vs IRG ---------------------------------------------------------------------------------
if (INPUT_ARGS$ISGvsIRG) {
	# This uses all differentially expressed genes
	# - Try to predict the direction of change (i.e. the sign)
	isgData <- FinalData %>% 
		filter(ExpressionCategory != 'NonDE')
	
	isgData$ExpressionCategory <- factor(isgData$ExpressionCategory, levels = c('Stimulated', 'Repressed')) # Make Stimulated the 'positive' value for logistic regression
	
	runid <- paste(c('ISGvsIRG', OutNameBase), collapse = '-')
	save_object(isgData, runid = runid, type = 'Data')
	
	# Training:
	train_model(isgData, runid = runid)
}



## ---- ISG vs other -------------------------------------------------------------------------------
if (INPUT_ARGS$ISGvsRandom) {
	# Compare ISG and NonDE only:
	isgData <- FinalData %>% 
		filter(ExpressionCategory != 'Repressed')
	
	isgData$ExpressionCategory <- factor(isgData$ExpressionCategory, levels = c('Stimulated', 'NonDE')) 
	
	runid <- paste(c('ISGvsRandom', OutNameBase), collapse = '-')
	save_object(isgData, runid = runid, type = 'Data')
	
	# Training:
	train_model(isgData, runid = runid)
}



## ---- IRG vs other -------------------------------------------------------------------------------
if (INPUT_ARGS$IRGvsRandom) {
	# Compare IRG and NonDE only:
	irgData <- FinalData %>% 
		filter(ExpressionCategory != 'Stimulated')
	
	irgData$ExpressionCategory <- factor(irgData$ExpressionCategory, levels = c('Repressed', 'NonDE')) 
	
	runid <- paste(c('IRGvsRandom', OutNameBase), collapse = '-')
	save_object(irgData, runid = runid, type = 'Data')
	
	# Training:
	train_model(irgData, runid = runid)
}



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Clean up -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
registerDoSEQ()
stopCluster(cl)

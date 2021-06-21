#!Rscript

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##   ISGpredict - Part 1
##     Merge data
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(rprojroot)
library(tidyverse)
library(argparse)

# Read user input:
parser <- ArgumentParser()
parser$add_argument('--host', nargs='+',
										help = 'Host to process. More precisely, if the data file takes the form "human_cdna_new_dat_fpkm_dups.txt", then "human" is the host')
INPUT_ARGS <- parser$parse_args()


# OTHER CONSTANTS:
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
part1Dir <- file.path(rootDir, 'Data', 'Part1_Dinucs')
part2Dir <- file.path(rootDir, 'Data', 'Part2_Codons')

source(UtilsFile)

# Read metadata:
metaDataPath <- file.path(rootDir, 'Data', 'ISG_identities.csv')
metaData <- read_csv(metaDataPath, na = c('NA', '-'))

# Read data for all included hosts:
part1Coding <- read_all_data(INPUT_ARGS$host, suffix = '_cds_new_dat_fpkm_dups.txt', directory = part1Dir)
part1Gene <- read_all_data(INPUT_ARGS$host, suffix = '_cdna_new_dat_fpkm_dups.txt', directory = part1Dir)

part2 <- read_all_data(INPUT_ARGS$host, suffix = '_cds_new_cpb_dat_dups.txt', directory = part2Dir)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Clean up metadata --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
sppMapping <- c('Bos taurus' = 'cow',
								'Canis familiaris' = 'dog',
								'Equus caballus' = 'horse',
								'Gallus gallus' = 'chicken',
								'Homo sapiens' = 'human',
								'Myotis lucifugus' = 'microbat',
								'Ovis aries' = 'sheep',
								'Pteropus vampyrus' = 'megabat',
								'Rattus norvegicus' = 'rat',
								'Sus scrofa' = 'pig')

expressionMapping <- c('up_regulated' = 'Stimulated',
											 'down_regulated' = 'Repressed',
											 'not_differentially_expressed' = 'NonDE')

metaData <- metaData %>% 
	select(HOST = Species, Gene = `ENSEMBL ID`, ExpressionCategory = Expression, 
				 Log2FC = `log2 Fold Change`, FDR = FDR) %>% 
	mutate(HOST = sppMapping[HOST],
				 ExpressionCategory = expressionMapping[ExpressionCategory]) %>% 
	distinct()
	
if (any(is.na(metaData$HOST))) stop('Mapping of host names failed. Check metadata.')
if (any(is.na(metaData$ExpressionCategory))) stop('Mapping of expression categories failed. Check metadata.')


if ('biomart' %in% INPUT_ARGS$host) {
	# Duplicate human data for joins to the biomart dataset:
	extraData <- metaData %>% 
		filter(HOST == 'human') %>% 
		mutate(HOST = 'biomart')
	
	metaData <- bind_rows(metaData, extraData)
}


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Rename columns to id data source -----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Columns to keep / rename
# - NOTE: currently we are dropping codon-pair bias columns (present in part2 data)
# - NOTE 2: Dropping transcript IDs - data has just one entry per gene, and transcript's used are
#           not always the same between gene and coding datasets

nameColumns <- c('^HOST$', '^Gene$') %>% 
	paste(collapse = '|')

part1Columns <- c(nameColumns, NUC_CONTENT, DINUC_BIAS, DINUC_PERCENT) %>% 
	paste(collapse = '|')

part2Columns <- c(nameColumns, BRIDGE_DINUC, NON_BRIDGE_DINUC, AMINO_BIAS, CODON_BIAS) %>% 
	paste(collapse = '|')


# All other columns renamed to identify source (coding vs entire gene):
#  - For part 2, we also have to process the SeqName column to extract Gene and Transcript IDs
part1Coding <- part1Coding %>% 
	select(matches(part1Columns, ignore.case = F)) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_Coding'))

part1Gene <- part1Gene %>% 
	select(matches(part1Columns, ignore.case = F)) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_EntireGene'))

part2 <- part2 %>% 
	mutate(SeqName = str_remove(SeqName, '^>')) %>% 
	separate(SeqName, into = c('Gene', 'Transcript'), sep = '\\|') %>% 
	select(matches(part2Columns, ignore.case = F)) %>% 
	rename_at(vars(-matches(nameColumns)), ~ paste0(., '_Coding'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Merge     ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
joinedData <- part1Coding %>% 
	full_join(part1Gene, by = c('HOST', 'Gene')) %>% 
	full_join(part2, by = c('HOST', 'Gene')) %>%
	left_join(metaData, by = c('HOST', 'Gene'))

# Check for duplication:
stopifnot(length(unique(joinedData$Gene)) == nrow(joinedData))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output     ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if (length(INPUT_ARGS$host) == 1) {
	outName <- paste0('InputData_', INPUT_ARGS$host, '.rds')
} else {
	warning('More than one host supplied - output will be named "InputData_AllHosts.rds"')  # These might not actually be all hosts in the analysis, but we have no way to check
	outName <- 'InputData_AllHosts.rds'
}

outPath <- file.path(rootDir, 'CalculatedData', outName)
saveRDS(joinedData, file = outPath)

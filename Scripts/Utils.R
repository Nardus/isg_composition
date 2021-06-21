## ISG prediction: Util functions


read_data <- function(host, suffix, directory) {
	## Read a single data file, fixing formatting issues if needed
	filename <- paste0(host, suffix) %>% 
		file.path(directory, .)
	
	# Part 2 files cause parsing issues because the first line does not end in a tab (while all 
	# other lines do)
	# - Check if file needs fixing and correct this before proceeding:
	rawLines <- readLines(filename)
	correctedLine <- paste0(rawLines[[1]], '\t')
	
	if (str_ends(rawLines[[2]], '\t') & !(str_ends(rawLines[[1]], '\t')))
		rawLines[[1]] <- correctedLine
	
	# Finally parse and return the data:
	read_tsv(rawLines)
}

read_all_data <- function(hosts, ...) {
	## Read data from multiple hosts, returning a merged data frame with a new column HOST 
	# identifying the host
	# ... are arguments to read_data, and must be named
	result <- lapply(hosts, read_data, ...)
	names(result) <- hosts
	
	bind_rows(result, .id = 'HOST')
}


read_rds_data <- function(host, prefix, folder = 'CalculatedData') {
	# Read RDS files by host, when all hosts are in the same folder
	dataFile <- paste0(prefix, host, '.rds')
	dataPath <- file.path(ROOT_DIR, folder, dataFile)
	
	readRDS(dataPath)
}


read_fits <- function(host, prefix, folder = 'Output') {
	# Read RDS files by host, when each host represents a subfolder
	dataFile <- paste0(prefix, host, '.rds')
	dataPath <- file.path(ROOT_DIR, folder, host, dataFile)
	
	readRDS(dataPath)
}


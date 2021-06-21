#
# Utility functions for plots
# 

HOST_COLOURS <- c(human = '#537AA7', rat = '#EC8D33', cow = '#DD5259', sheep = '#78B6B2',
									pig = '#5D9C4A', horse = '#EFCA4E', dog = '#B67FA4', megabat = '#F1A1AE',
									microbat = '#A37A64', chicken = '#C3B8B3')

save_plot <- function(p, basename, basedir = ROOT_DIR, plotdir = file.path('Output', 'Plots'), ...) {
	# Save plots as pdf
	# - p can be a list, in which case multiple plots are saved to the same file
	# - otherwise, p is assumed to a valid single plot (e.g. from ggplot)
	outpath <- file.path(basedir, plotdir)
	
	if (!dir.exists(outpath)) {
		dir.create(outpath)
	}
	
	pdf(file.path(outpath, paste0(basename, '.pdf')), ...)
	if ('list' %in% class(p)) {
		lapply(p, print)
		
	} else {
		print(p)
	}
	dev.off()
}


## ------ For reading run data -------------------------------------------------------------------------
read_data <- function(host, suffix, directory) {
	filename <- paste0(host, suffix) %>% 
		file.path(directory, .)
	
	read_csv(filename)
}


read_fits <- function(host, prefix, dataDir, hostSubDir = FALSE) {
	# READ FITS:
	# - if hostSubDir is true, look for a subdirectory inside the main data dir for
	#   each host
	
	filename <- paste0(prefix, host, '.rds')
	
	
	if (hostSubDir) {
		filename <- file.path(ROOT_DIR, dataDir, host, filename)
	} else {
		filename <- file.path(ROOT_DIR, dataDir, filename)
	}
	
	readRDS(filename)
}


read_by_runid <- function(host, runid, type = 'Data', folder = '') {
	prefix <- paste(type, runid, '', sep = '_') 
	dataDir <- file.path('Output', folder, runid)
	
	read_fits(host, prefix, dataDir)
}


## ------ Other utils for processing feature names: ------------------------------------------------
escape_special_chars <- function(x) {
	# Escape names containing special chars by adding backticks
	# Useful for recreating the column names used by readr:: funtions
	
	allowedX <- gsub('[._]', '', x)  # Remove dots and underscores, which are allowed (should not be escaped)
	
	startsSpecial <- grepl('^[_[:digit:]]', x)  # Can't start with an underscore or a letter
	containsSpecial <- grepl('[[:punct:]]', allowedX)
	
	# Escape and return:
	ifelse(startsSpecial | containsSpecial,
				 paste0('`', x, '`'),
				 x)
}


add_feature_class <- function(data) {
	# Add a FeatureClass column to data
	# (expects data to a have a column named 'Feature')
	
	data %>% 
		mutate(FeatureClass = gsub('_[[:alnum:]]+$', '', Feature),  # Remove trailing '_Coding'/'_EntireGene'/etc.
					 FeatureClass = gsub('(^br)|(^NonBr)', '', FeatureClass), # Strip of leading bridge / non-bridge identifiers
					 FeatureClass = gsub('^([[:upper:]])([[:upper:]])%$', '\\1p\\2', FeatureClass))  # Convert names like 'AG%' to 'ApG'
}



## ------ Variable importance: ---------------------------------------------------------------------
get_top_N <- function(varimps, Ntop) {
	# Get top features by Rank
	top <- varimps %>%
		filter(Rank <= Ntop) %>%
		arrange(Rank)
	
	top$Feature <- factor(top$Feature, levels = rev(unique(top$Feature)))
	
	return(top)
}


plot_top_vars <- function(varimps, yval, xlabel, ylabel = 'Feature', Ntop = 50) {
	# Plot the top N variables, with their importance. 
	yval <- enquo(yval)
	
	top50 <- get_top_N(varimps, Ntop)
	
	ggplot(top50, aes(x = Feature, y = !!yval, fill = IsCpG)) +
		geom_bar(stat = 'identity', position = 'dodge') +
		scale_fill_manual(guide = F, values = c('TRUE' = 'red', 'FALSE' = 'grey50')) +
		coord_flip() +
		labs(x = ylabel, y = xlabel) +  # swapped because of coord_flip
		PLOT_THEME
}


# For multiclass predictions, we have multiple importance values:
plot_top_vars_multiclass <- function(varimps, yval, xlabel, ylabel = 'Feature', Ntop = 50) {
	# Plot the top N variables, with their importance
	yval <- enquo(yval)
	
	top50 <- get_top_N(varimps, Ntop)
	
	ggplot(top50, aes(x = reorder_within(Feature, by = !!yval, within = Label), y = !!yval, fill = IsCpG)) +
		geom_bar(stat = 'identity', position = 'dodge') +
		scale_fill_manual(guide = F, values = c('TRUE' = 'red', 'FALSE' = 'grey50')) +
		facet_wrap(~ Label, scales = 'free', nrow = 1) +
		scale_x_reordered() +
		coord_flip() +
		labs(x = ylabel, y = xlabel) +  # swapped because of coord_flip
		PLOT_THEME
}


add_title <- function(p, host, baseTitle) {
	fullTitle <- paste(baseTitle, host, sep = ': ')
	
	p + 
		ggtitle(fullTitle)
}


plot_CpG_ranks <- function(fit, data, plotTitle = NULL) {
	# Plot the ranks of CpG features in different hosts:
	# - Ranks plotted for individual genes
	
	# Get gene-level varimps:
	varimps <- get_shap_contributions(fit, collapse = TRUE) %>% 
		left_join(data, by = c('RowID' = 'rowIndex')) %>% 
		group_by(Host, Gene) %>% 
		mutate(Rank = rank(-abs(SHAP), ties.method = 'min')) %>%   # Largest Shap value ranked 1st
		filter(FeatureClass == 'CpG') %>% 
		mutate(Rank = ifelse(SHAP == 0, Inf, Rank))  # Unused features will have importance of 0, so their rank becomes arbitrary
	
	# Set y-axis limits to include 1 (and not 0):
	get_lims <- function(lims) {
		res <- pretty(lims)
		res <- res[res > 1]
		c(1, res)
	}
	
	# Plot:
	p <- ggplot(varimps, aes(x = Host, y = Rank, colour = ExpressionCategory)) +
		geom_boxplot(fill = NA) +
		scale_colour_brewer(palette = 'Set1') +
		scale_y_reverse(expand = c(0.02, 0), breaks = get_lims) +
		labs(title = plotTitle, colour = 'Label') +
		#facet_grid(. ~ Host, scales = 'free_x', space = 'free_x') +
		PLOT_THEME +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
					strip.text = element_blank(),
					panel.spacing.x = unit(0, 'lines'))
	
	# Return data and plot:
	varimps <- varimps %>% 
		select(RowID, Host, Gene, ExpressionCategory, SHAP, Rank)
	
	list(cpg_ranks = varimps, plot = p)
}


get_shap_contributions <- function(trainedModel, collapse = FALSE, get_interactions = FALSE) {
	# Extract SHAP scores for individual vars and observations
	# - Based on code in xgb.plot.shap()
	x <- as.matrix(trainedModel$trainingData[, -1])
	colnames(x) <- escape_special_chars(colnames(x))  # predict.xgb checks that column names are identical to model$feature_names
	
	model <- trainedModel$finalModel
	
	shap <- xgboost:::predict.xgb.Booster(model, x, predcontrib = TRUE, approxcontrib = FALSE, predinteraction = FALSE)
	
	# For multiclass models:
	# - `shap` will contain log odds of being in each class (versus in the other categories)
	if (is.list(shap)) {
		names(shap) <- levels(trainedModel$trainingData$.outcome)
		
		shap <- lapply(shap, data.frame, check.names = FALSE) %>% 
			bind_rows(.id = 'Label')
		
	} else {
		shap <- data.frame(shap, check.names = FALSE)
		shap$Label <- levels(trainedModel$trainingData$.outcome)[2]  # 1st level is the baseline, so the shap values are log odds of being label 2
	}
	
	# Clean up:
	shap <- shap %>% 
		mutate(RowID = 1:n()) %>% 
		select(-BIAS) %>% 
		gather(-RowID, -Label, key = 'Feature', value = 'SHAP') %>% 
		mutate(Feature = gsub('`', '', Feature))
	
	# Collapse to a single value per gene for each feature class (if requested)
	#   - Note that while SHAP values are additive (by definition), this is 
	#   merely an assumption about the structure of the model
	if (collapse) {
		shap <- shap %>% 
			add_feature_class() %>% 
			group_by(RowID, FeatureClass, Label) %>% 
			summarise(SHAP = sum(SHAP))
	}
	
	# Return:
	return(shap)
}


get_varimp_shap <- function(trainedModel, collapse = FALSE, keepDirection = FALSE) {
	shapContribs <- get_shap_contributions(trainedModel, collapse = collapse)
	
	if (collapse) {
		# Rename so all other functions will work regardless of whether features have been collapsed or not
		shapContribs <- shapContribs %>% 
			rename(Feature = FeatureClass)
	}
	
	if (keepDirection) {
		transf <- function(x) x  # keep as is
	} else {
		transf <- abs   # get absolute value
	}
	
	varImpShap <- shapContribs %>% 
		group_by(Feature, Label) %>%     # Group by Label: only applies for multiclass models (otherwise label is identical for all genes)
		summarise(MeanSHAP = mean(transf(SHAP))) %>% 
		ungroup() %>% 
		add_feature_class() %>% 
		mutate(IsCpG = FeatureClass == 'CpG',
					 Rank = rank(-MeanSHAP, ties.method = 'min')) %>%   # Highes Shap value ranked 1st
		arrange(desc(MeanSHAP))
	
	varImpShap$Feature <- factor(varImpShap$Feature, levels = unique(varImpShap$Feature))
	
	varImpShap
}





#
plot_overall_importance <- function(varimps, xlabel = 'Relative importance', ylabel = 'Feature', 
																		title = NULL, Ntop = 50, normalise = TRUE, rescale = FALSE) {
	# Plot importance accross all hosts
	# - If importance estimates come from separate models, set normalise = T (default) to ensure 
	#   importance estimates are comparable
	# - Setting rescale = TRUE will re-scale importance estimates to sum to 1, allowing easier
	#   comparison between plots.
	
	# Rank by overall importance to get the most important features in general
	if (normalise) {
		varimps <- varimps %>% 
			group_by(Host) %>% 
			mutate(RelativeImportance = MeanSHAP / max(MeanSHAP))  # Relative importance, to make values comparable across spp
	} else {
		varimps <- varimps %>% 
			mutate(RelativeImportance = MeanSHAP)
	}
	
	
	totalImportance <- varimps %>% 
		group_by(Feature) %>% 
		summarise(TotalImportance = sum(RelativeImportance)) %>% 
		ungroup() %>% 
		mutate(Rank = rank(-TotalImportance, ties.method = 'min'))
	
	topFeatures <- varimps %>% 
		select(-Rank) %>% 
		left_join(totalImportance, by = 'Feature') %>% 
		get_top_N(Ntop)
	
	if (rescale) {
		topFeatures <- topFeatures %>% 
			ungroup() %>% 
			mutate(RelativeImportance = RelativeImportance / sum(RelativeImportance))
	}
	
	print(sum(topFeatures$RelativeImportance))
	
	mainPlot <- ggplot(topFeatures, aes(x = Feature, y = RelativeImportance, fill = Host)) +
		geom_col(position = 'stack', colour = 'grey20') +
		coord_flip() +
		scale_fill_brewer(palette = 'Spectral', drop = FALSE) +
		labs(title = title, x = ylabel, y = xlabel) +  # swapped because of coord_flip
		PLOT_THEME
	
	# Also plot species missing in each bar:
	missingHosts <- topFeatures %>% 
		filter(RelativeImportance == 0) %>% 
		distinct(Feature, Host) %>% 
		mutate(RelativeImportance = 1)
	
	if (nrow(missingHosts) == 0) {
		# Nothing to add
		return(mainPlot)
	}
	
	mainPlot <- mainPlot + 
		guides(fill = FALSE)
	
	sidePlot <- ggplot(missingHosts, aes(x = Feature, y = RelativeImportance, fill = Host)) +
		geom_col(position = 'stack', colour = 'grey20', size = 0.1) +
		scale_fill_brewer(palette = 'Spectral', drop = FALSE) +
		scale_y_continuous(name = 'Not present:', position = 'right') +
		coord_flip() +
		PLOT_THEME +
		theme(axis.text = element_blank(),
					axis.title.y = element_blank(),
					axis.ticks.x = element_blank(),
					panel.grid = element_blank())
	
	
	ggarrange(mainPlot, sidePlot, widths = c(5, 1), draw = FALSE)
}

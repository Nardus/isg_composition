#
# Plots for runs containing multiple hosts
#


library(rprojroot)
library(tidyverse)
library(xgboost)
library(caret)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(zoo)
library(tidytext)
library(ClustOfVar)


# Constants:
ROOT_DIR <- find_rstudio_root_file()

PLOT_THEME <- theme_bw() +
	theme(panel.grid = element_blank(),
				panel.border = element_rect(size = 1))

CLASS_COLOURS <- brewer.pal(3, 'Dark2')
names(CLASS_COLOURS) <- c('Stimulated', 'NonDE', 'Repressed')


# Utility functions:
UtilsFile <- file.path(ROOT_DIR, 'Scripts', 'Utils.R') # Contains utils for reading raw data
source(UtilsFile)

PlotUtilsFile <- file.path(ROOT_DIR, 'Scripts', 'PlotUtils.R')
source(PlotUtilsFile)

save_plot2 <- function(...) {
	save_plot(..., plotdir = file.path('Output', 'AllHostsCombined', 'Plots'))
}

HOSTS <- names(HOST_COLOURS)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Input data
inputData.ISGvsIRG <- read_by_runid('AllHosts', runid = 'ISGvsIRG-Top50', type = 'Data', folder = 'AllHostsCombined')
inputData.StimulatedVsRandom <- read_by_runid('AllHosts', runid = 'ISGvsRandom-Top50', type = 'Data', folder = 'AllHostsCombined')
inputData.RepressedVsRandom <- read_by_runid('AllHosts', runid = 'IRGvsRandom-Top50', type = 'Data', folder = 'AllHostsCombined')

inputData.ISGvsIRG <- inputData.ISGvsIRG %>% 
	mutate(rowIndex = 1:n(),
				 Host = factor(HOST, levels = HOSTS))

inputData.StimulatedVsRandom <- inputData.StimulatedVsRandom %>% 
	mutate(rowIndex = 1:n(),
				 Host = factor(HOST, levels = HOSTS))

inputData.RepressedVsRandom <- inputData.RepressedVsRandom %>% 
	mutate(rowIndex = 1:n(),
				 Host = factor(HOST, levels = HOSTS))

# Model fits
fit.ISGvsIRG <- read_by_runid('AllHosts', runid = 'ISGvsIRG-Top50', type = 'Fit', folder = 'AllHostsCombined')
fit.StimulatedVsRandom <- read_by_runid('AllHosts', runid = 'ISGvsRandom-Top50', type = 'Fit', folder = 'AllHostsCombined')
fit.RepressedVsRandom <- read_by_runid('AllHosts', runid = 'IRGvsRandom-Top50', type = 'Fit', folder = 'AllHostsCombined')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_accuracy <- function(fit, data, title = '', by_host = FALSE) {
	# Calculate accuracy for each host in each iteration:
	# Variance across iterations used to guage instability (by showing error bars)
	preds <- fit$pred %>% 
		left_join(data, by = 'rowIndex') %>% 
		separate(Resample, into = c('Fold', 'Replicate'), sep = '\\.') %>% 
		ungroup()
	
	if (by_host) {
		accuracy <- preds %>% 
			group_by(obs, Replicate, Host) %>% 
			summarise(Accuracy = sum(pred == obs) / n()) %>% 
			group_by(Host, obs) %>% 
			summarise(MeanAccuracy = mean(Accuracy),   # Per class accuracy, so really sensitivity and specificity
								LwrAccuracy = quantile(Accuracy, probs = 0.05/2),
								UprAccuracy = quantile(Accuracy, probs = 1 - 0.05/2))
		
	} else {
		accuracy <- preds %>% 
			group_by(obs, Replicate) %>% 
			summarise(Accuracy = sum(pred == obs) / n()) %>% 
			group_by(obs) %>% 
			summarise(MeanAccuracy = mean(Accuracy),
								LwrAccuracy = quantile(Accuracy, probs = 0.05/2),
								UprAccuracy = quantile(Accuracy, probs = 1 - 0.05/2))
	}
	
	# Plot:
	if (by_host) {
		p <- ggplot(accuracy, aes(x = Host, y = MeanAccuracy, fill = obs)) +
			xlab('Species')
	} else {
		p <- ggplot(accuracy, aes(x = obs, y = MeanAccuracy, fill = obs)) +
			xlab('Class')
	}
	
	p <- p +
		geom_bar(stat = 'identity', position = 'dodge') +
		geom_errorbar(aes(ymin = LwrAccuracy, ymax = UprAccuracy), position = 'dodge', colour = 'grey20') +
		geom_hline(yintercept = 0.5, linetype = 2, colour = 'grey30') +
		labs(title = title, y = 'Accuracy', fill = 'Class') +
		ylim(0, 1) +
		scale_fill_manual(values = CLASS_COLOURS) +
		guides(fill = FALSE) +
		PLOT_THEME +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
	
	list(p = p, data = accuracy)
}


make_class_legend <- function(colours) {
	# Colours should be a named vector
	dt <- data.frame(x = 1, y = 1,
									 col = names(colours))
	
	p <- ggplot(dt, aes(x = x, y = y, fill = col)) +
		geom_col() +
		scale_fill_manual(values = CLASS_COLOURS, name = 'Class')
	
	get_legend(p)
}

p1 <- plot_accuracy(fit.ISGvsIRG, inputData.ISGvsIRG)

p2 <- plot_accuracy(fit.StimulatedVsRandom, inputData.StimulatedVsRandom)

p3 <- plot_accuracy(fit.RepressedVsRandom, inputData.RepressedVsRandom)


# Combine:
combined <- plot_grid(p1$p, p2$p, p3$p,
											nrow = 1,
											labels = c('A', 'B', 'C'), hjust = 0.2,
											align = 'v', axis = 'b') +
	theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))


save_plot2(combined, 'OverallAccuracy_Combined', height = 3, width = 7.5)


# Output plot data:
list(A = p1$data, B = p2$data, C = p3$data) %>% 
	bind_rows(.id = 'Panel') %>% 
	write_excel_csv(file.path('Output', 'AllHostsCombined', 'Plots', 'OverallAccuracy_Combined.csv'))


raw_preds.ISGvsIRG <- fit.ISGvsIRG$pred %>% 
	left_join(inputData.ISGvsIRG, by = 'rowIndex') %>% 
	separate(Resample, into = c('Fold', 'Replicate'), sep = '\\.') %>% 
	select(Gene, Species = Host, Observed = obs, Fold, Replicate, Prediction = pred) %>% 
	arrange(Fold, Replicate, Species, Gene)

raw_preds.StimulatedVsRandom <- fit.StimulatedVsRandom$pred %>% 
	left_join(inputData.StimulatedVsRandom, by = 'rowIndex') %>% 
	separate(Resample, into = c('Fold', 'Replicate'), sep = '\\.') %>% 
	select(Gene, Species = Host, Observed = obs, Fold, Replicate, Prediction = pred) %>% 
	arrange(Fold, Replicate, Species, Gene)

raw_preds.RepressedVsRandom <- fit.RepressedVsRandom$pred %>% 
	left_join(inputData.RepressedVsRandom, by = 'rowIndex') %>% 
	separate(Resample, into = c('Fold', 'Replicate'), sep = '\\.') %>% 
	select(Gene, Species = Host, Observed = obs, Fold, Replicate, Prediction = pred) %>% 
	arrange(Fold, Replicate, Species, Gene)

write_excel_csv(raw_preds.ISGvsIRG, file.path('Output', 'AllHostsCombined', 'Plots', 'OverallAccuracy_RawPredictions_ISGvsIRG.csv'))
write_excel_csv(raw_preds.StimulatedVsRandom, file.path('Output', 'AllHostsCombined', 'Plots', 'OverallAccuracy_RawPredictions_StimulatedVsRandom.csv'))
write_excel_csv(raw_preds.RepressedVsRandom, file.path('Output', 'AllHostsCombined', 'Plots', 'OverallAccuracy_RawPredictions_RepressedVsRandom.csv'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy by gene ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_accuracy_genelevel <- function(fit, data, title = '') {
	# Get accuracy by gene:
	accuracy <- fit$pred %>% 
		left_join(data, by = 'rowIndex') %>% 
		group_by(rowIndex, Host, obs) %>% 
		summarise(Accuracy = sum(pred == obs) / n()) %>%  # Per class accuracy, so really sensitivity and specificity
		ungroup() %>% 
		mutate(OutlierHigh = Accuracy > quantile(Accuracy, .75) + 1.50*IQR(Accuracy),
					 OutlierLow = Accuracy < quantile(Accuracy, .25) - 1.50*IQR(Accuracy),
					 OutlierShape = ifelse(OutlierHigh | OutlierLow, 'A', NA))
	
	ggplot(accuracy, aes(x = obs, y = Accuracy)) +
		geom_boxplot(outlier.shape = NA) +
		geom_jitter(aes(shape = OutlierShape), height = 0.01, width = 0.1) +  # Show oultiers, but jitter them
		geom_hline(yintercept = 0.5, linetype = 2) +
		labs(title = title, x = NULL, y = 'Mean accuracy') +
		guides(shape = FALSE) +
		ylim(0, 1) +
		PLOT_THEME +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


p1 <- plot_accuracy_genelevel(fit.ISGvsIRG, inputData.ISGvsIRG, 
															title = 'ISG versus IRG: overall accuracy (all hosts in a single run)')

p2 <- plot_accuracy_genelevel(fit.StimulatedVsRandom, inputData.StimulatedVsRandom, 
															title = 'ISG versus Random: overall accuracy (all hosts in a single run)')

p3 <- plot_accuracy_genelevel(fit.RepressedVsRandom, inputData.RepressedVsRandom, 
															title = 'IRG versus Random: overall accuracy (all hosts in a single run)')



save_plot2(p1, 'Accuracy_GeneLevel_ISGvsIRG')
save_plot2(p2, 'Accuracy_GeneLevel_ISGvsRandom')
save_plot2(p3, 'Accuracy_GeneLevel_IRGvsRandom')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy across hosts ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_accuracy_byhost <- function(fit, data, title = '') {
	# Get accuracy by gene:
	accuracy <- fit$pred %>% 
		left_join(data, by = 'rowIndex') %>% 
		group_by(Host, obs) %>% 
		summarise(Accuracy = sum(pred == obs) / n()) %>%  # Per class accuracy, so really sensitivity and specificity
		ungroup() %>% 
		mutate(OutlierHigh = Accuracy > quantile(Accuracy, .75) + 1.50*IQR(Accuracy),
					 OutlierLow = Accuracy < quantile(Accuracy, .25) - 1.50*IQR(Accuracy),
					 OutlierShape = ifelse(OutlierHigh | OutlierLow, 'A', NA))
	
	ggplot(accuracy, aes(x = obs, y = Accuracy)) +
		geom_boxplot(outlier.shape = NA) +
		geom_jitter(aes(colour = Host), height = 0, width = 0.05) +
		geom_hline(yintercept = 0.5, linetype = 2) +
		labs(title = title, x = NULL, y = 'Mean accuracy') +
		guides(shape = FALSE) +
		ylim(0, 1) +
		PLOT_THEME +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


p1 <- plot_accuracy_byhost(fit.ISGvsIRG, inputData.ISGvsIRG, 
													 title = 'ISG versus IRG: accuracy by host (all hosts in a single run)')

p2 <- plot_accuracy_byhost(fit.StimulatedVsRandom, inputData.StimulatedVsRandom, 
													 title = 'ISG versus Random: accuracy by host (all hosts in a single run)')

p3 <- plot_accuracy_byhost(fit.RepressedVsRandom, inputData.RepressedVsRandom, 
													 title = 'IRG versus Random: accuracy by host (all hosts in a single run)')


save_plot2(p1, 'Accuracy_ByHost_ISGvsIRG')
save_plot2(p2, 'Accuracy_ByHost_ISGvsRandom')
save_plot2(p3, 'Accuracy_ByHost_IRGvsRandom')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Variable importance ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_varimp_shap_multihost <- function(trainedModel, inputData, collapse = FALSE, keepDirection = FALSE) {
	shapContribs <- get_shap_contributions(trainedModel, collapse = collapse) %>% 
		left_join(inputData, by = c('RowID' = 'rowIndex'))
	
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
	
	# Get the mean Shapley value for each feature class across all genes
	# - For display, this value is divided among the reservoirs based on the effect size (shapley value) of 
	#   the feature in their genes
	# - TODO: Legend will need careful wording since this may be confusing
	varImpShap <- shapContribs %>% 
		group_by(Feature, Label) %>%        # Group by Label: only applies for multiclass models (otherwise label is identical for all genes)
		mutate(OverallMeanSHAP = mean(transf(SHAP))) %>% 
		group_by(Host, Feature, Label) %>%     
		summarise(HostMeanSHAP = mean(transf(SHAP)),
							OverallMeanSHAP = unique(OverallMeanSHAP)) %>% 
		ungroup() %>% 
		add_feature_class() %>% 
		mutate(IsCpG = FeatureClass == 'CpG',
					 HostContribution = OverallMeanSHAP * (HostMeanSHAP/OverallMeanSHAP),  # Proportion of bar which should be give over to this host
					 HostContribution = HostContribution / sum(HostContribution, na.rm = T))     # Re-scale, so total across all features is 1
	
	# Add ranks:
	ranks <- varImpShap %>% 
		distinct(Feature, OverallMeanSHAP) %>% 
		mutate(Rank = rank(-OverallMeanSHAP, ties.method = 'min')) %>%   # Highest Shap value ranked 1st)
		select(-OverallMeanSHAP)
	
	varImpShap <- varImpShap %>% 
		left_join(ranks, by = 'Feature') %>% 
		arrange(OverallMeanSHAP)
	
	varImpShap$Feature <- factor(varImpShap$Feature, levels = unique(varImpShap$Feature))
	
	varImpShap
}


plot_importance_by_host <- function(varimps, Ntop, xlabel = 'Relative importance', ylabel = 'Feature', 
																	 title = NULL) {
	varimps %>% 
		filter(Rank <= Ntop) %>% 
		ggplot(aes(x = Feature, y = HostContribution, fill = Host)) +
			geom_col(position = 'stack') +
			coord_flip() +
			#scale_fill_manual(values = HOST_COLOURS) +
			scale_fill_brewer(palette = 'Spectral') +
			labs(title = title, x = ylabel, y = xlabel) +  # swapped because of coord_flip
			PLOT_THEME
}


# Plots:
varImpShap.ISGvsIRG <- get_varimp_shap_multihost(fit.ISGvsIRG, inputData.ISGvsIRG, collapse = TRUE)
varImpShap.StimulatedVsRandom <- get_varimp_shap_multihost(fit.StimulatedVsRandom, inputData.StimulatedVsRandom, collapse = TRUE)
varImpShap.RepressedVsRandom <- get_varimp_shap_multihost(fit.RepressedVsRandom, inputData.RepressedVsRandom, collapse = TRUE)

	

p1 <- plot_importance_by_host(varImpShap.ISGvsIRG, Ntop = 15, xlabel = 'Relative importance') +
	guides(fill = F) +
	ylim(0, 0.2)

p2 <- plot_importance_by_host(varImpShap.StimulatedVsRandom, Ntop = 15, xlabel = 'Relative importance') +
	guides(fill = F) +
	xlab(NULL) +
	ylim(0, 0.2)

p3 <- plot_importance_by_host(varImpShap.RepressedVsRandom, Ntop = 15, xlabel = 'Relative importance') +
	xlab(NULL) +
	theme(legend.box.margin = margin(b = 30)) +
	ylim(0, 0.2)


# Combine:
shared_legend <- get_legend(p3)
p3 <- p3 + guides(fill = F)

p_combined <- plot_grid(p1, p2, p3, shared_legend,
					rel_widths = c(1.1, 1, 1, 0.5),  # First plot contains shared y-axis label
					labels = c('D', 'E', 'F', ''), vjust = 0.5, hjust = 0.2,
					nrow = 1, align = 'v', axis = 'b') +
	theme(plot.margin = margin(10.5, 5.5, 5.5, 5.5))


save_plot2(p_combined, 'ImportanceAcrossSpp_Combined', height = 3.5, width = 7.5)

#save_plot2(p1, 'ImportanceAcrossSpp_ISGvsIRG_stacked', height = 3.5, width = 7.5/3.5)
#save_plot2(p2, 'ImportanceAcrossSpp_ISGvsRandom_stacked', height = 3.5, width = 7.5/3.5)
#save_plot2(p3, 'ImportanceAcrossSpp_IRGvsRandom_stacked', height = 3.5, width = 7.5/3.5)
#save_plot2(plot_grid(shared_legend), 'ImportanceAcrossSpp_Legend', height = 3.5, width = 7.5/3.5*0.5)


# Save raw data from plots
importance.ISGvsIRG <- varImpShap.ISGvsIRG %>% 
	mutate(HostContribution = if_else(is.na(HostContribution), 0, HostContribution)) %>% 
	arrange(-HostContribution) %>% 
	select(Species = Host, Feature, FeatureImportance_AllSpecies = OverallMeanSHAP, FeatureImportance_ThisSpecies = HostMeanSHAP, RelativeImportance_ThisSpecies = HostContribution)

importance.StimulatedVsRandom <- varImpShap.StimulatedVsRandom %>% 
	mutate(HostContribution = if_else(is.na(HostContribution), 0, HostContribution)) %>% 
	arrange(-HostContribution) %>% 
	select(Species = Host, Feature, FeatureImportance_AllSpecies = OverallMeanSHAP, FeatureImportance_ThisSpecies = HostMeanSHAP, RelativeImportance_ThisSpecies = HostContribution)

importance.RepressedVsRandom <- varImpShap.RepressedVsRandom %>% 
	mutate(HostContribution = if_else(is.na(HostContribution), 0, HostContribution)) %>% 
	arrange(-HostContribution) %>% 
	select(Species = Host, Feature, FeatureImportance_AllSpecies = OverallMeanSHAP, FeatureImportance_ThisSpecies = HostMeanSHAP, RelativeImportance_ThisSpecies = HostContribution)


write_excel_csv(importance.ISGvsIRG, file.path('Output', 'AllHostsCombined', 'Plots', 'Importance_ISGvsIRG.csv'))
write_excel_csv(importance.StimulatedVsRandom, file.path('Output', 'AllHostsCombined', 'Plots', 'Importance_StimulatedVsRandom.csv'))
write_excel_csv(importance.RepressedVsRandom, file.path('Output', 'AllHostsCombined', 'Plots', 'Importance_RepressedVsRandom.csv'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- CpG rank by gene ---------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cpg_rank.ISGvsIRG <- plot_CpG_ranks(fit.ISGvsIRG, inputData.ISGvsIRG, plotTitle = NULL)
cpg_rank.StimulatedVsRandom <- plot_CpG_ranks(fit.StimulatedVsRandom, inputData.StimulatedVsRandom, plotTitle = NULL)
cpg_rank.RepressedVsRandom <- plot_CpG_ranks(fit.RepressedVsRandom, inputData.RepressedVsRandom, plotTitle = NULL)

p1 <- cpg_rank.ISGvsIRG$plot +
	scale_colour_manual(values = CLASS_COLOURS, guide = FALSE)

p2 <- cpg_rank.StimulatedVsRandom$plot +
	scale_colour_manual(values = CLASS_COLOURS, guide = FALSE)

p3 <- cpg_rank.RepressedVsRandom$plot +
	scale_colour_manual(values = CLASS_COLOURS, guide = FALSE)


# Combine:
shared_legend <- make_class_legend(CLASS_COLOURS)

combined <- plot_grid(p1, p2, p3, shared_legend,
					nrow = 1, rel_widths = c(1, 1, 1, 0.5),
					labels = c('A', 'B', 'C'),
					align = 'v', axis = 't')

save_plot2(combined, 'CpGranks_Combined', height = 5)

#save_plot2(p1, 'CpGranks_ISGvsIRG')
#save_plot2(p2, 'CpGranks_ISGvsRandom')
#save_plot2(p3, 'CpGranks_IRGvsRandom')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Correlation of top features with CpG -------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_top_feature_correlations <- function(data, varimps, ntop = 50, title = 'Correlation with CpG_EntireGene') {
	topFeatures <- varimps %>% 
		distinct(Feature, IsCpG, OverallMeanSHAP, Rank) %>% 
		filter(Rank <= ntop)
	
	cors <- data %>% 
		gather(-Host, -HOST, -Gene, -Log2FC, -FDR, -ExpressionCategory, -rowIndex,
					 key = 'Feature', value = 'Value') %>% 
		right_join(topFeatures, by = 'Feature') %>% 
		group_by(Feature, IsCpG, Rank, OverallMeanSHAP) %>% 
		summarise(Cor_with_CpG = cor(Value, data$CpG_EntireGene, method = 'spearman', use = 'pairwise.complete.obs')) %>% 
		ungroup() %>% 
		arrange(-Rank) %>% 
		mutate(Feature = factor(Feature, levels = unique(Feature)))
	
	fplot <- ggplot(cors, aes(x = Feature, y = OverallMeanSHAP, fill = IsCpG)) +
		geom_col() +
		coord_flip() +
		scale_fill_brewer(palette = 'Set1', direction = -1, guide = FALSE) +
		ggtitle(title)
	
	cplot <- ggplot(cors, aes(x = Feature, y = Cor_with_CpG, fill = IsCpG)) +
		geom_col() +
		geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, alpha = 0.5) +
		ylim(c(-1, 1)) +
		coord_flip() +
		scale_fill_brewer(palette = 'Set1', direction = -1, guide = FALSE) +
		theme(axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					plot.background = element_blank())
	
	plot_grid(fplot, cplot, nrow = 1, align = 'h')
}

# Need uncollapsed versions of varimps:
varImpShap.ISGvsIRG <- get_varimp_shap_multihost(fit.ISGvsIRG, inputData.ISGvsIRG, collapse = FALSE)
varImpShap.StimulatedVsRandom <- get_varimp_shap_multihost(fit.StimulatedVsRandom, inputData.StimulatedVsRandom, collapse = FALSE)
varImpShap.RepressedVsRandom <- get_varimp_shap_multihost(fit.RepressedVsRandom, inputData.RepressedVsRandom, collapse = FALSE)

p1 <- plot_top_feature_correlations(data = inputData.ISGvsIRG,
																		varimps = varImpShap.ISGvsIRG,
																		title = 'ISG vs IRG (correlation with CpG_EntireGene)')

p2 <- plot_top_feature_correlations(data = inputData.StimulatedVsRandom,
																		varimps = varImpShap.StimulatedVsRandom,
																		title = 'ISG vs Random (correlation with CpG_EntireGene)')

p3 <- plot_top_feature_correlations(data = inputData.RepressedVsRandom,
																		varimps = varImpShap.RepressedVsRandom,
																		title = 'IRG vs Random (correlation with CpG_EntireGene)')

save_plot2(p1, 'Correlations_ISGvsIRG')
save_plot2(p2, 'Correlations_ISGvsRandom')
save_plot2(p3, 'Correlations_IRGvsRandom')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Direction for individual effects (partial effect-like plots) -------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
plot_feature_effects <- function(trainedModel, inputData, varimps, ntop = 12, ylabel = 'Log odds', title = NULL) {
	# Get feature values
	inputData <- inputData %>% 
		gather(-HOST, -Host, -Gene, -Log2FC, -FDR, -ExpressionCategory, -rowIndex,
					 key = 'Feature', value = 'Value')
	
	# Effects
	topFeatures <- varimps %>% 
		distinct(Feature, Rank) %>% 
		filter(Rank <= ntop) %>% 
		arrange(Rank) %>% 
		.$Feature
	
	shapContribs <- get_shap_contributions(trainedModel, collapse = FALSE) %>% 
		left_join(inputData, by = c('RowID' = 'rowIndex', 'Feature' = 'Feature')) %>% 
		filter(Feature %in% topFeatures) %>% 
		mutate(Feature = factor(Feature, levels = topFeatures))
	
	# Plot
	ggplot(shapContribs, aes(x = Value, y = SHAP)) +
		geom_point(aes(colour = ExpressionCategory), size = 0.1, shape = 1) +
		#geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) +
		geom_hline(yintercept = 0, linetype = 2) +
		facet_wrap(Feature ~ ., scales = 'free_x') +
		scale_colour_brewer(palette = 'Set1', direction = -1) +
		labs(x = 'Feature value', y = ylabel, title = title) +
		PLOT_THEME
}

p1 <- plot_feature_effects(trainedModel = fit.ISGvsIRG, 
													 inputData = inputData.ISGvsIRG, 
													 varimps = varImpShap.ISGvsIRG,
													 ylabel = 'Log odds ISG', title = 'ISG vs IRG', ntop = 25)

p2 <- plot_feature_effects(trainedModel = fit.StimulatedVsRandom, 
													 inputData = inputData.StimulatedVsRandom, 
													 varimps = varImpShap.StimulatedVsRandom,
													 ylabel = 'Log odds ISG', title = 'ISG vs Random', ntop = 25)

p3 <- plot_feature_effects(trainedModel = fit.RepressedVsRandom, 
													 inputData = inputData.RepressedVsRandom, 
													 varimps = varImpShap.RepressedVsRandom,
													 ylabel = 'Log odds IRG', title = 'IRG vs Random', ntop = 25)


save_plot2(p1, 'EffectDirections_ISGvsIRG', width = 11, height = 11)
save_plot2(p2, 'EffectDirections_ISGvsRandom', width = 11, height = 11)
save_plot2(p3, 'EffectDirections_IRGvsRandom', width = 11, height = 11)



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Variable importance, but also show direction -----------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
varImpShap.ISGvsIRG <- get_varimp_shap_multihost(fit.ISGvsIRG, 
																								 inputData.ISGvsIRG, 
																								 collapse = TRUE, 
																								 keepDirection = TRUE)

varImpShap.StimulatedVsRandom <- get_varimp_shap_multihost(fit.StimulatedVsRandom, 
																													 inputData.StimulatedVsRandom, 
																													 collapse = TRUE, 
																													 keepDirection = TRUE)

varImpShap.RepressedVsRandom <- get_varimp_shap_multihost(fit.RepressedVsRandom, 
																													inputData.RepressedVsRandom, 
																													collapse = TRUE,
																													keepDirection = TRUE)


p1 <- plot_overall_importance(varImpShap.ISGvsIRG, Ntop = 30, normalise = FALSE,
															xlabel = 'Importance (total effect on log odds)',
															title = 'Importance across species: ISG vs IRG\n(all hosts in a single run)')

p2 <- plot_overall_importance(varImpShap.StimulatedVsRandom, Ntop = 30, normalise = FALSE,
															xlabel = 'Importance (total effect on log odds)',
															title = 'Importance across species: ISG vs Random\n(all hosts in a single run)')

p3 <- plot_overall_importance(varImpShap.RepressedVsRandom, Ntop = 30, normalise = FALSE,
															xlabel = 'Importance (total effect on log odds)',
															title = 'Importance across species: IRG vs Random\n(all hosts in a single run)')


save_plot2(p1, 'Importance_ISGvsIRG')
save_plot2(p2, 'Importance_ISGvsRandom')
save_plot2(p3, 'Importance_IRGvsRandom')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- CpG raw data -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
all_raw_data <- inputData.StimulatedVsRandom %>% 
	filter(ExpressionCategory == 'NonDE') %>% 
	bind_rows(inputData.ISGvsIRG) %>% 
	mutate(ExpressionCategory = factor(ExpressionCategory, levels = c('Stimulated', 'NonDE', 'Repressed')))

p1 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = CpG_EntireGene, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	ggtitle('All hosts combined')

p2 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = CpG_EntireGene, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	facet_wrap(~ Host) +
	ggtitle('By host')

p3 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = `S-Bias_Coding`, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	ggtitle('All hosts combined')

p4 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = `S-Bias_Coding`, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	facet_wrap(~ Host) +
	ggtitle('By host')

p5 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = `TGC-Bias_Coding`, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	ggtitle('All hosts combined')

p6 <- ggplot(all_raw_data, aes(x = ExpressionCategory, y = `TGC-Bias_Coding`, colour = ExpressionCategory)) +
	geom_violin() +
	geom_boxplot(width = 0.3) +
	scale_color_brewer(palette = 'Set1', guide = FALSE) +
	facet_wrap(~ Host) +
	ggtitle('By host')

save_plot2(p1, 'RawData_CpG_Combined')
save_plot2(p2, 'RawData_CpG_ByHost')

save_plot2(p3, 'RawData_S-Bias_Combined')
save_plot2(p4, 'RawData_S-Bias_ByHost')

save_plot2(p5, 'RawData_TGC-Bias_Combined')
save_plot2(p6, 'RawData_TGC-Bias_ByHost')
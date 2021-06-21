#!Rscript

## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
##   ISGpredict - Part 3
##     Plot variable importances
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(rprojroot)
library(tidyverse)
library(xgboost)
library(caret)
library(egg)
library(zoo)
library(tidytext)
library(ClustOfVar)
library(argparse)


# Read user input:
parser <- ArgumentParser()
parser$add_argument('host', type = 'character', nargs = '+',
										help = 'one or more hosts to include in the plots. Order of hosts determines plotting order')

INPUT_ARGS <- parser$parse_args()  #c('human', 'rat', 'megabat', 'microbat', 'cow', 'pig', 'sheep', 'horse', 'dog', 'chicken'))
HOSTS <- INPUT_ARGS$host

ROOT_DIR <- find_rstudio_root_file()

PLOT_THEME <- theme_bw()


# Utility functions:
UtilsFile <- file.path(ROOT_DIR, 'Scripts', 'Utils.R') # Contains utils for reading raw data
source(UtilsFile)

PlotUtilsFile <- file.path(ROOT_DIR, 'Scripts', 'PlotUtils.R')
source(PlotUtilsFile)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Input data
inputData.ISGvsIRG <- lapply(HOSTS, read_by_runid, runid = 'ISGvsIRG-Top50', type = 'Data')
inputData.StimulatedVsRandom <- lapply('human', read_by_runid, runid = 'ISGvsRandom-Top50', type = 'Data')
inputData.RepressedVsRandom <- lapply('human', read_by_runid, runid = 'IRGvsRandom-Top50', type = 'Data')


names(inputData.ISGvsIRG) <- HOSTS
names(inputData.StimulatedVsRandom) <- 'human'
names(inputData.RepressedVsRandom) <- 'human'


# Model fits
fit.ISGvsIRG <- lapply(HOSTS, read_by_runid, runid = 'ISGvsIRG-Top50', type = 'Fit')
fit.StimulatedVsRandom <- lapply('human', read_by_runid, runid = 'ISGvsRandom-Top50', type = 'Fit')
fit.RepressedVsRandom <- lapply('human', read_by_runid, runid = 'IRGvsRandom-Top50', type = 'Fit')

names(fit.ISGvsIRG) <- HOSTS
names(fit.StimulatedVsRandom) <- 'human'
names(fit.RepressedVsRandom) <- 'human'



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Expected predictiveness of different CpG vars: ---------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
cpgDists.vsIRG <- inputData.ISGvsIRG %>% 
	bind_rows(.id = 'Host') %>% 
	select(Host, ExpressionCategory, contains('CpG'), contains('CG%')) %>% 
	gather(-Host, -ExpressionCategory, key = 'Feature', value = 'FeatureValue') %>% 
	mutate(RowID = 1:n()) %>% 
	spread(key = ExpressionCategory, value = FeatureValue) %>% 
	group_by(Host, Feature) %>% 
	do(Dist = as.vector(dist(cbind(.$Stimulated, .$Repressed)))) %>% 
	unnest()
	

ggplot(cpgDists.vsIRG, aes(x = Feature, y = Dist, colour = Host)) +
	geom_boxplot() +
	coord_flip()





# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy -----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_accuracy <- function(trainedModel) {
	finalPreds <- trainedModel$pred
	
	finalPreds %>% 
		group_by(Resample, obs) %>% 
		summarise(Accuracy = sum(pred == obs) / n()) %>%  # Per class accuracy, so really sensitivity and specificity
		group_by(obs) %>% 
		summarise(MeanAccuracy = mean(Accuracy))
}

accuracy.ISGvsIRG <- lapply(fit.ISGvsIRG, get_accuracy) %>% 
	bind_rows(.id = 'Host')

accuracy.StimulatedVsRandom <- lapply(fit.StimulatedVsRandom, get_accuracy) %>% 
	bind_rows(.id = 'Host')

accuracy.RepressedVsRandom <- lapply(fit.RepressedVsRandom, get_accuracy) %>% 
	bind_rows(.id = 'Host')

# Force plotting order to match input order:
accuracy.ISGvsIRG$Host <- factor(accuracy.ISGvsIRG$Host, levels = HOSTS)
accuracy.StimulatedVsRandom$Host <- factor(accuracy.StimulatedVsRandom$Host, levels = HOSTS)
accuracy.RepressedVsRandom$Host <- factor(accuracy.RepressedVsRandom$Host, levels = HOSTS)


# Plots:
p1 <- ggplot(accuracy.ISGvsIRG, aes(x = Host, y = MeanAccuracy)) +
	geom_bar(stat = 'identity') +
	geom_hline(yintercept = 0.5, linetype = 2) +
	facet_wrap(~ obs) +
	labs(title = 'Top 50 ISG versus IRG', x = NULL, y = 'Mean accuracy') +
	ylim(0, 1) +
	PLOT_THEME +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p2 <- ggplot(accuracy.StimulatedVsRandom, aes(x = Host, y = MeanAccuracy)) +
	geom_bar(stat = 'identity') +
	geom_hline(yintercept = 0.5, linetype = 2) +
	facet_wrap(~ obs) +
	labs(title = 'Top 50 ISG versus random genes', x = NULL, y = 'Mean accuracy') +
	ylim(0, 1) +
	PLOT_THEME +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p3 <- ggplot(accuracy.RepressedVsRandom, aes(x = Host, y = MeanAccuracy)) +
	geom_bar(stat = 'identity') +
	geom_hline(yintercept = 0.5, linetype = 2) +
	facet_wrap(~ obs) +
	labs(title = 'Top 50 IRG versus random genes', x = NULL, y = 'Mean accuracy') +
	ylim(0, 1) +
	PLOT_THEME +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))



save_plot(p1, 'OverallAccuracy_ISGvsIRG')
save_plot(p2, 'OverallAccuracy_ISGvsRandom')
save_plot(p3, 'OverallAccuracy_IRGvsRandom')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Accuracy relative to expression level ------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_gene_level_accuracy <- function(host, modelList, inputDataList, metadataList = TranscriptData) {
	# Get gene-level accuracy, and add original raw data plus metadata about the transcripts
	trainedModel <- modelList[[host]]
	inputData <- inputDataList[[host]]
	
	accuracyData <- trainedModel$pred %>%
		group_by(rowIndex) %>%
		summarise(Accuracy = sum(pred == obs) / n()) %>%
		ungroup()
	
	inputData %>%
		ungroup() %>%
		mutate(rowIndex = 1:n()) %>%
		left_join(accuracyData, by = 'rowIndex')
}


plot_accuracy_by_expression <- function(modelList, inputDataList, windowSize = 100) {
	accuracyData <- lapply(names(modelList), get_gene_level_accuracy, 
												 modelList = modelList, inputDataList = inputDataList) %>% 
		bind_rows() %>% 
		rename(Host = HOST) %>% 
		group_by(ExpressionCategory, Host) %>% 
		mutate(Rank = rank(abs(Log2FC), ties.method = 'min')) %>% 
		ungroup()
	
	# Calculate accuracy over a sliding window, as a way to show the raw data
	accuracySliding <- accuracyData %>% 
		arrange(Rank) %>% 
		group_by(Host, ExpressionCategory) %>%
		mutate(Accuracy = rollapply(Accuracy, FUN = mean, width = windowSize, partial = TRUE)) %>% 
		filter(Rank <= (max(Rank) - windowSize))  # Trim end, where window has increasingly fewer data
	
	
	# Fix plotting order:
	accuracyData$Host <- factor(accuracyData$Host, levels = HOSTS)
	
	# Plot:
	ggplot(accuracyData, aes(x = Rank, y = Accuracy, colour = Host, fill = Host)) +
		geom_line(data = accuracySliding, alpha = 0.5) +
		geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"), se = TRUE) +
		geom_hline(yintercept = 0.5, linetype = 2) +
		facet_grid(Host ~ ExpressionCategory, scales = 'free') +
		labs(x = 'Gene rank', y = 'Mean accuracy') +
		PLOT_THEME
}




p1 <- plot_accuracy_by_expression(modelList = fit.ISGvsIRG, inputDataList = inputData.ISGvsIRG, windowSize = 10) +
	ggtitle('ISG vs IRG: Genes ranked by magnitude of fold change')
	

p2 <- plot_accuracy_by_expression(modelList = fit.StimulatedVsRandom, inputDataList = inputData.StimulatedVsRandom, windowSize = 10) +
	ggtitle('Top 50 ISG vs Random')


p3 <- plot_accuracy_by_expression(modelList = fit.RepressedVsRandom, inputDataList = inputData.RepressedVsRandom, windowSize = 10) +
	ggtitle('Top 50 IRG vs Random')




save_plot(p1, 'AccuracyByTranscripType_ISGvsIRG')
save_plot(p2, 'AccuracyByTranscriptType_ISGvsRandom')
save_plot(p3, 'AccuracyByTranscriptType_IRGvsRandom')




# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Variable importance ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Calculations:
varImpShap.ISGvsIRG <- lapply(fit.ISGvsIRG, get_varimp_shap, collapse = TRUE)
varImpShap.StimulatedVsRandom <- lapply(fit.StimulatedVsRandom, get_varimp_shap, collapse = TRUE)
varImpShap.RepressedVsRandom <- lapply(fit.RepressedVsRandom, get_varimp_shap, collapse = TRUE)


# Plots:
p1List <- lapply(varImpShap.ISGvsIRG, plot_top_vars, yval = MeanSHAP, 
								 xlabel = 'Change in log odds', ylabel = 'Feature class')
p1List <- mapply(add_title, p = p1List, host = names(p1List),
								 MoreArgs = list(baseTitle = 'Variable importance (ISG vs IRG)'), SIMPLIFY = F)

p2List <- lapply(varImpShap.StimulatedVsRandom, plot_top_vars, yval = MeanSHAP, 
								 xlabel = 'Change in log odds', ylabel = 'Feature class')
p2List <- mapply(add_title, p = p2List, host = names(p2List),
								 MoreArgs = list(baseTitle = 'Variable importance (ISG vs Random)'), SIMPLIFY = F)

p3List <- lapply(varImpShap.RepressedVsRandom, plot_top_vars, yval = MeanSHAP, 
								 xlabel = 'Change in log odds', ylabel = 'Feature class')
p3List <- mapply(add_title, p = p2List, host = names(p2List),
								 MoreArgs = list(baseTitle = 'Variable importance (IRG vs Random)'), SIMPLIFY = F)



p4 <- varImpShap.ISGvsIRG %>% 
	bind_rows(.id = 'Host') %>% 
	plot_CpG_ranks(rankedby = MeanSHAP, plotTitle = 'Relative position of CpG (ISG vs IRG)') +
	ylim(0, 10)  # Make sure Inf is clearly distinguishable from other ranks

p5 <- varImpShap.StimulatedVsRandom %>% 
	bind_rows(.id = 'Host') %>% 
	plot_CpG_ranks(rankedby = MeanSHAP, plotTitle = 'Relative position of CpG (ISG vs Random)')

p6 <- varImpShap.RepressedVsRandom %>% 
	bind_rows(.id = 'Host') %>% 
	plot_CpG_ranks(rankedby = MeanSHAP, plotTitle = 'Relative position of CpG (IRG vs Random)')


save_plot(p1List, 'VarImp-SHAP_Top50_ISGvsIRG', height = 15)
save_plot(p2List, 'VarImp-SHAP_Top50_ISGvsRandom', height = 15)
save_plot(p3List, 'VarImp-SHAP_Top50_IRGvsISG', height = 15)

save_plot(p4, 'VarImp-SHAP_CpGRanks_ISGvsIRG')
save_plot(p5, 'VarImp-SHAP_CpGRanks_ISGvsRandom')
save_plot(p6, 'VarImp-SHAP_CpGRanks_IRGvsRandom')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Variable importance across all species -----------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
overalImp <- varImpShap.ISGvsIRG %>% 
	bind_rows(.id = 'Host')

overalImp$Host <- factor(overalImp$Host, levels = HOSTS)


stackedPlot <- plot_overall_importance(overalImp, Ntop = 30, xlabel = 'Importance',
																			 title = 'Importance across species: ISG vs IRG (top 50)\n')


save_plot(stackedPlot, 'ImportanceAcrossSpp_ISGvsIRG_stacked')



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Does CpG importance decline by distance? ---------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
get_varimp_raw <- function(trainedModel, collapse = FALSE) {
	# Like get_varimp_shap, but summarise at sample level
	# TODO: can be combined with above by adding a summarise option?
	
	shapContribs <- get_shap_contributions(trainedModel, collapse = collapse)
	
	if (collapse) {
		# Rename so all other functions will work regardless of whether features have been collapsed or not
		shapContribs <- shapContribs %>% 
			rename(Feature = FeatureClass)
	}
	
	varImpShap <- shapContribs %>% 
		add_feature_class() %>% 
		mutate(IsCpG = FeatureClass == 'CpG') %>% 
		group_by(RowID) %>%   # Features ranked for each sample
		mutate(Rank = rank(-abs(SHAP), ties.method = 'min')) %>%   # Rank largest SHAP values 1st
		arrange(desc(abs(SHAP)))
	
	varImpShap$Feature <- factor(varImpShap$Feature, levels = unique(varImpShap$Feature))
	
	varImpShap
}

plot_CpG_rank_genelevel <- function(varimps, plotTitle = NULL, addFacet = TRUE) {
	# Plot the ranks of CpG features as a function of gene expression
	varimps <- filter(varimps, IsCpG) %>% 
		group_by(ExpressionCategory) %>% 
		mutate(Top50 = if_else(ExpressionCategory == 'Stimulated', GeneRank >= (max(GeneRank)-50), 
													 if_else(ExpressionCategory == 'Repressed', GeneRank <= 50,
													 				 FALSE)),
					 Rank = ifelse(SHAP == 0, Inf, Rank))  # Unused vars do not have a defined rank...
	
	p <- ggplot(varimps, aes(x = Log2FC, y = Rank)) +
		geom_point(aes(colour = Top50), alpha = 0.5) +
		#geom_rug(aes(colour = Top50), sides = 'l', position = 'jitter') +
		#geom_smooth(method = 'gam', formula = y ~ s(x, bs = "tp")) +
		scale_colour_manual(guide = F, values = c('TRUE' = 'red', 'FALSE' = 'grey40')) +
		scale_fill_manual(guide = F, values = c('TRUE' = 'red', 'FALSE' = 'grey40')) +
		scale_y_reverse() +
		labs(x = 'Log2FC', y = 'CpG rank') +
		ggtitle(plotTitle) +
		PLOT_THEME +
		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
	
	if (addFacet) {
		p <- p + facet_grid(. ~ ExpressionCategory, scales = 'free_x', space = 'free')
	}
	
	p
}


make_genelevel_plot <- function(host, runid, plotTitle, facetOrder = NULL, addFacet = TRUE) {
	inputData <- read_by_runid(host, runid)
	fit <- read_by_runid(host, runid, type = 'Fit')
	
	inputData <- inputData %>% 
		mutate(RowID = 1:n()) %>% 
		group_by(ExpressionCategory) %>% 
		mutate(GeneRank = rank(Log2FC, ties.method = 'min'))  # Can't rank by absolute fold change here (does not make sense for random genes)
	
	varImpSampleLevel <- get_varimp_raw(fit, collapse = TRUE) %>% 
		left_join(inputData, by = 'RowID')
	
	# Allow custom ordering of facets (e.g. to broadly match x-axis order):
	if (!is.null(facetOrder)) {
		varImpSampleLevel$ExpressionCategory <- factor(as.character(varImpSampleLevel$ExpressionCategory), levels = facetOrder)
	}
	
	plot_CpG_rank_genelevel(varImpSampleLevel, plotTitle, addFacet = addFacet)
}


# Use all genes / human here:
p1 <- make_genelevel_plot('human', 'ISGvsRandom-AllGenes', plotTitle = 'Relative position of CpG (ISG vs Random - human, all genes)',
													facetOrder = c('NonDE', 'Stimulated'))

p2 <- make_genelevel_plot('human', 'ISGvsIRG-AllGenes', plotTitle = 'Relative position of CpG (ISG vs IRG - human, all genes)', addFacet = FALSE)

p3 <- make_genelevel_plot('human', 'IRGvsRandom-AllGenes', plotTitle = 'Relative position of CpG (IRG vs Random - human, all genes)', 
													facetOrder = c('Repressed', 'NonDE'))

save_plot(p1, 'CpG_Importance_by_FoldChange-ISGvsRandom')
save_plot(p2, 'CpG_Importance_by_FoldChange-ISGvsIRG')
save_plot(p3, 'CpG_Importance_by_FoldChange-IRGvsRandom')

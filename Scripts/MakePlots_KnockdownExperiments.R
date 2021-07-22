# 
# Plot predictions for A549 knockdown experiments
# 

library(argparse)

parser <- ArgumentParser()
parser$add_argument("model", type = "character",
										help = "base name of model to to use (e.g. 'ISGvsIRG')")
parser$add_argument("positive_class", type = "character",
										help = "label of positive class (e.g. 'ISG')")

INPUT_ARGS <- parser$parse_args()

# Other libraries
library(rprojroot)
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(ggplot2)
library(cowplot)

ROOT_DIR <- find_rstudio_root_file()
source(file.path(ROOT_DIR, 'Scripts', 'PlotUtils.R'))

set.seed(2342147)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Trained model
full_model_name <- sprintf("%s-Top50", INPUT_ARGS$model)
model_file_name <- sprintf("Fit_%s_AllHosts.rds", full_model_name)

trained_model <- file.path(ROOT_DIR, 'Output', 'AllHostsCombined', full_model_name, model_file_name) %>% 
	readRDS()

# Labels depend on which model is used:
probability_label <- sprintf("Predicted probability of being an %s", INPUT_ARGS$positive_class)
cpg_effect_label <- sprintf("CpG effect (on %s probability)", INPUT_ARGS$positive_class)
output_name <- sprintf("Output/AllHostsCombined/Plots/cpg+probability_%s_model.pdf", INPUT_ARGS$model)
data_folder <- sprintf("Output/predictions/%s", INPUT_ARGS$model)
dir.create(data_folder, recursive = TRUE)


# IRGs from fibroblast experiment (Shaw et al. 2017)
top_irgs <- file.path('Data', 'SelectedGenes', 'human_Top50_Formatted.csv') %>% 
	read_csv(col_types = cols(.default = 'c')) %>% 
	select(EnsemblID = .data$`Ensembl gene ID`, ExpressionCategory = .data$Class)

# IRGs consistent across all experiments:
consistent_irgs <- file.path("Data", "Other_Experiments", "consistent_irgs.csv") %>% 
	read_csv(col_types = cols(.default = "c"))

# Random genes NOT used during training:
holdout_random <- file.path("Data", "Other_Experiments", "Random100_nonDE_CPM1_A549_typeI_4h.csv") %>% 
	read_csv(col_types = cols(.default = "c")) %>% 
	select(.data$EnsemblID) %>% 
	mutate(Category = "Holdout Random") %>% 
	sample_n(size = 50)


# Data for all human genes
features <- file.path(ROOT_DIR, 'CalculatedData', 'InputData_biomart.rds') %>% 
	readRDS()

all_gene_metadata <- read.csv("Data/ISG_identities.csv") %>% 
	filter(.data$Species == "Homo sapiens") %>% 
	distinct(EnsemblID = .data$ENSEMBL.ID, Gene)

top_irgs <- top_irgs %>% 
	left_join(all_gene_metadata, by = "EnsemblID")

holdout_random <- holdout_random %>% 
	left_join(all_gene_metadata, by = "EnsemblID")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- General predictions ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
colnames(features) <- escape_special_chars(colnames(features))  # predict.xgb checks that column names are identical to model$feature_names

feature_mat <- features %>% 
	select(one_of(trained_model$finalModel$feature_names)) %>% 
	as.matrix()

preds <- data.frame(EnsemblID = features$Gene,
										predicted_prob = xgboost:::predict.xgb.Booster(trained_model$finalModel, feature_mat))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- CpG effect ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Effect of CpG features on predictions for all genes:
shap <- xgboost:::predict.xgb.Booster(trained_model$finalModel, feature_mat, predcontrib = TRUE, approxcontrib = FALSE, predinteraction = FALSE)

shap <- data.frame(shap, check.names = FALSE) %>% 
	mutate(EnsemblID = features$Gene,
				 Log2FC = features$Log2FC,
				 ExpressionCategory = features$ExpressionCategory)


# Summarise to higher level "feature class" (i.e. all CpG measures combined)
shap_long <- shap %>% 
	mutate(RowID = 1:n()) %>% 
	select(-.data$BIAS) %>% 
	gather(-.data$RowID, -.data$EnsemblID, -.data$Log2FC, -.data$ExpressionCategory, key = 'Feature', value = 'SHAP') %>% 
	mutate(Feature = gsub('`', '', .data$Feature)) %>% 
	add_feature_class() 


shap_summary <- shap_long %>% 
	filter(.data$FeatureClass == "CpG") %>% 
	group_by(.data$RowID, .data$FeatureClass, .data$EnsemblID, .data$Log2FC, .data$ExpressionCategory) %>% 
	summarise(SHAP = sum(.data$SHAP)) %>% 
	ungroup()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Overview plots: top 50 IRGs ----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
stopifnot(!any(holdout_random$EnsemblID %in% top_irgs$EnsemblID))  # Should be no overlap with existing categories

top_irgs <- top_irgs %>% 
	mutate(Category = case_when(.data$ExpressionCategory == "IRG" ~ "Top 50 IRGs",
															.data$ExpressionCategory == "ISG" ~ "Top 50 ISGs",
															.data$ExpressionCategory == "Random" ~ "Random")) %>% 
	bind_rows(holdout_random) %>% 
	left_join(shap_summary, by = "EnsemblID") %>% 
	left_join(preds, by = "EnsemblID")


overview_prob <- ggplot(top_irgs, aes(x = Category, y = predicted_prob)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = probability_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

overview_cpg <- ggplot(top_irgs, aes(x = Category, y = SHAP)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = cpg_effect_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Consistent IRGs across all 3 RNA-seq experiments (with KO response) ------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
stopifnot(all(consistent_irgs$ensemble_id %in% shap_summary$EnsemblID))
stopifnot(all(consistent_irgs$ensemble_id %in% preds$EnsemblID))

ko_consistent <- consistent_irgs %>% 
	select(Gene = .data$gene,
				 EnsemblID = .data$ensemble_id,
				 Category = .data$ko_experiment) %>% 
	left_join(shap_summary, by = "EnsemblID") %>% 
	left_join(preds, by = "EnsemblID") %>% 
	bind_rows(top_irgs)


koc_prob <- ggplot(ko_consistent, aes(x = Category, y = predicted_prob)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = probability_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

koc_cpg <- ggplot(ko_consistent, aes(x = Category, y = SHAP)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = cpg_effect_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Save predictions
ko_consistent %>% 
	select(-.data$FeatureClass, -.data$Log2FC, -starts_with("ExpressionCategory"), -.data$RowID,
				 cpg_effect = .data$SHAP) %>% 
	write_excel_csv(file.path(data_folder, "5C.csv"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Knockout experiment ------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# random / IRGs (no guide) / IRGs (KO)
# TODO: unclear where random genes should come from (can we really assume the original set are all still random here?)
koexp_control_data <- read_excel("Data/Other_Experiments/EdgeRCas9ClonesIFNome_withFDR_forPrism.xlsx",
																 sheet = "EdgeRCas9ClonesIFNome_withFDR_p") %>% 
	mutate(Treatment = "no guide") %>% 
	select(.data$EnsemblID, .data$Treatment, label = .data$Rank)

koexp_ko_data <- read_excel("Data/Other_Experiments/EdgeR_KOClonesIFNome_withFDR_forPrisim.xlsx",
														sheet = "EdgeR_KOClonesIFNome_withFDR") %>% 
	mutate(Treatment = "KO") %>% 
	filter(!is.na(.data$EnsemblID) & .data$EnsemblID != "NA") %>% 
	select(.data$EnsemblID, .data$Treatment, label = .data$Rank)


ko_experiment <- bind_rows(koexp_control_data, koexp_ko_data) %>% 
	mutate(Category = case_when(.data$label == "a-Top50-ISGs" ~ "ISG",
															.data$label == "c-Top50-IRGs" ~ "IRG"),
				 Label = sprintf("%s (%s)", .data$Category, .data$Treatment))

koexp_random <- holdout_random %>% 
	mutate(Label = .data$Category)


ko_experiment <- ko_experiment %>%
	left_join(all_gene_metadata, by = "EnsemblID") %>% 
	bind_rows(koexp_random) %>% 
	left_join(shap_summary, by = "EnsemblID") %>% 
	left_join(preds, by = "EnsemblID")


koexp_prob <- ggplot(ko_experiment, aes(x = Label, y = predicted_prob)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = probability_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

koexp_cpg <- ggplot(ko_experiment, aes(x = Label, y = SHAP)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = cpg_effect_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Save predictions
ko_experiment %>% 
	select(-.data$FeatureClass, -.data$ExpressionCategory, -.data$Log2FC, -.data$RowID, -.data$label,
				 cpg_effect = .data$SHAP) %>% 
	write_excel_csv(file.path(data_folder, "4D_new.csv"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Exogenous ZAP - all genes -----------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
exo_zap_genes <- read.csv("Data/Other_Experiments/Fig_5I_Top50-Up-Down-Random.csv") %>% 
	select(Gene = .data$Gene_Name, .data$EnsemblID, Category = .data$Rank) %>% 
	bind_rows(holdout_random) %>% 
	mutate(Category = factor(.data$Category, levels = c("TOP50 UP", "Random50", "Holdout Random", "TOP50 DOWN"))) %>% 
	left_join(shap_summary, by = "EnsemblID") %>% 
	left_join(preds, by = "EnsemblID")


exozap_prob <- ggplot(exo_zap_genes, aes(x = Category, y = predicted_prob)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = probability_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

exozap_cpg <- ggplot(exo_zap_genes, aes(x = Category, y = SHAP)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = cpg_effect_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Save predictions
exo_zap_genes %>% 
	select(-.data$FeatureClass, -.data$ExpressionCategory, -.data$Log2FC, -.data$RowID,
				 cpg_effect = .data$SHAP) %>% 
	write_excel_csv(file.path(data_folder, "5F.csv"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Overexpression - consistent IRGs -----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
exozap_consistent <- consistent_irgs %>% 
	select(Gene = .data$gene, 
				 EnsemblID = .data$ensemble_id, 
				 Category = .data$exogenous_zap_experiment) %>% 
	filter(!is.na(.data$Category)) %>% 
	left_join(shap_summary, by = "EnsemblID") %>% 
	left_join(preds, by = "EnsemblID") %>% 
	bind_rows(top_irgs)


exoc_prob <- ggplot(exozap_consistent, aes(x = Category, y = predicted_prob)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = probability_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

exoc_cpg <- ggplot(exozap_consistent, aes(x = Category, y = SHAP)) +
	geom_boxplot() +
	geom_jitter(height = 0, width = 0.1) +
	labs(x = NULL, y = cpg_effect_label) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# Save predictions
exozap_consistent %>% 
	select(-.data$FeatureClass, -starts_with("ExpressionCategory"), -.data$Log2FC, -.data$RowID,
				 cpg_effect = .data$SHAP) %>% 
	write_excel_csv(file.path(data_folder, "5H.csv"))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Combine plots and save ---------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
add_title <- function(p, t) {
	p + 
		ggtitle(t)
}

titles <- c("Fibroblast experiment", 
						"KO - consistent IRGs (5C)", 
						"KO - IRGs (4D)",
						"Exogenous ZAP (5F)",
						"Exogenous ZAP (5H)")

p_prob <- list(overview_prob, koc_prob, koexp_prob, exozap_prob, exoc_prob) %>% 
	mapply(add_title, p = ., t = titles, SIMPLIFY = FALSE) %>% 
	plot_grid(plotlist = ., nrow = 1)

p_cpg <- list(overview_cpg, koc_cpg, koexp_cpg, exozap_cpg, exoc_cpg) %>% 
	mapply(add_title, p = ., t = titles, SIMPLIFY = FALSE) %>% 
	plot_grid(plotlist = ., nrow = 1)


pdf(output_name, width = 15, height = 4)
	p_cpg
	p_prob
dev.off()

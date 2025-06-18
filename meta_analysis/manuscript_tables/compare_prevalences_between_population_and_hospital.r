library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

source("../meta_analysis_utils.r")

# Hospital based biobanks
# bbj, biome, mgbb, pmbb.

# Population level
# uk-biobank, all-of-us, genes-and-health, egcut.

# Mixed
# ccpm, gel

# I want to compare the hospital-based biobanks to the population-based ones,
# using Wilcoxon tests

# First, let's determine the case prevalences of each of the traits within each 
# of the biobanks

files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks",
	recursive=TRUE)
dt_info <- rbindlist(lapply(files, extract_file_info), fill=TRUE)
dt_info <- dt_info %>% mutate(dataset = gsub(".*/", "", dataset)) %>% 
	mutate(sampling_strategy = ifelse(
		dataset %in% c("bbj", "mgbb", "pmbb", "biome"),
		"Hospital/Health center–based", "Population-based"))
dt_info <- dt_info %>%
	mutate(sampling_strategy = ifelse(
		dataset %in% c("ccpm", "gel"), "Mixed", sampling_strategy))
dt_info <- dt_info %>% filter(binary) %>% mutate(
	Prevalence = as.integer(n_cases)/
	(as.integer(n_cases) + as.integer(n_controls))) %>%
mutate(
	Biobank = unlist(renaming_plot_biobank_list[dataset]),
	Description = unlist(renaming_phenotype_list[phenotype]))

fwrite(dt_info,
	"Tables/table_for_creation_of_prevalence_plots.tsv",
	sep="\t", quote=FALSE)

dt_info <- fread("Tables/table_for_creation_of_prevalence_plots.tsv") %>%
	filter(phenotype %in% 
		c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
		"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
		"CRF", "CRP", "EFRMB", "FemInf","Gout", "HDLC", "Height", "HF",
		"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
		"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
		"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
		"WHRBMI", "AlcCons", "ALT", "HipRep"))

dt_info <- dt_info %>% mutate(hospital = case_when(
	sampling_strategy == "Hospital/Health center–based" ~ TRUE,
	sampling_strategy == "Population-based" ~ FALSE,
	.default = NA)
)
dt_info_test <- dt_info %>% filter(phenotype != "HipRep")
# Then, perform the test, split by phenotype
results <- dt_info_test %>% group_by(phenotype) %>% 
	summarise(P_value = wilcox.test(
		Prevalence ~ hospital,
		paired=FALSE,
		alternative="less")$`p.value`)
print(results %>% arrange(P_value), n=50)

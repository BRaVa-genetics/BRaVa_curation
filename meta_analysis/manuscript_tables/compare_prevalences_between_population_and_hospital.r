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
	full.names=TRUE, recursive=TRUE)
files <- basename(grep("cleaned/gene", files, value=TRUE))
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
	Description = unlist(renaming_phenotype_list[phenotype])) %>%
	filter(phenotype %in% 
		c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
		"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
		"CRF", "CRP", "EFRMB", "FemInf","Gout", "HDLC", "Height", "HF",
		"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
		"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
		"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
		"WHRBMI", "AlcCons", "ALT", "HipRep"))

fwrite(dt_info,
	"Tables/table_for_creation_of_prevalence_plots.tsv",
	sep="\t", quote=FALSE)

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

# Create a cleaned table for the paper
dt_info <- fread("Tables/table_for_creation_of_prevalence_plots.tsv") %>%
	filter(phenotype %in% 
		c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
		"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
		"CRF", "CRP", "EFRMB", "FemInf","Gout", "HDLC", "Height", "HF",
		"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
		"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
		"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
		"WHRBMI", "AlcCons", "ALT", "HipRep"))

dt_info <- data.table(dt_info %>% 
	rename(`Phenotype ID` = phenotype, `Biobank ID` = dataset,
		`N cases` = n_cases, `N controls` = n_controls) %>%
	select(Description, `Phenotype ID`, Biobank, `Biobank ID`, `N cases`, `N controls`) %>%
	group_by(Biobank, Description, `Phenotype ID`, `Biobank ID`) %>% 
	summarise(
		`N cases` = sum(`N cases`),
		`N controls` = sum(`N controls`),
		Prevalence = sum(`N cases`) / 
		sum(`N cases` + `N controls`)))

setkeyv(dt_info, c("Description", "Biobank"))
fwrite(dt_info, "Tables/supplementary_table_9.tsv")

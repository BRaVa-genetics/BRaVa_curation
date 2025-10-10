#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(latex2exp)

# Loop over the data files that contain the case and control counts
# # Local
# files <- dir(path="../data/meta_analysis/gcloud", full.names=TRUE, recursive=TRUE)
# BMRC
files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs", full.names=TRUE, recursive=TRUE)
files <- gsub("^.*/", "", grep("/cleaned/", files, value=TRUE))
source("../meta_analysis_utils.r")

file_info <- rbindlist(lapply(files, extract_file_info), fill=TRUE)
file_info <- file_info %>% filter(phenotype %in% 
	c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
	"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
	"CRF", "CRP", "EFRMB", "FemInf","Gout", "HDLC", "Height", "HF",
	"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
	"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
	"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
	"WHRBMI", "AlcCons", "ALT", "HipRep"))
file_info <- file_info %>% filter(type == "gene")

# Now, split by phenotype as well and include the counts
phenotype_ancestry_information <- file_info %>% group_by(dataset, ancestry, phenotype) %>%
	summarise(
		n=as.integer(n),
		n_cases=as.integer(n_cases),
		n_controls=as.integer(n_controls), na.rm=TRUE
	)

phenotype_ancestry_information <- phenotype_ancestry_information %>%
	rename(
		`Biobank ID` = dataset,
		`Phenotype ID` = phenotype,
		Ancestry = ancestry,
		N = n,
		`N cases` = n_cases,
		`N controls` = n_controls) %>%
	mutate(
		Biobank = unlist(renaming_plot_biobank_list[`Biobank ID`]),
		Description = unlist(renaming_phenotype_list[`Phenotype ID`]),
		)

# Binary traits
phenotype_ancestry_information_binary <- phenotype_ancestry_information %>% filter(!is.na(`N cases`)) %>%
	select(Description, `Phenotype ID`, Ancestry, Biobank, `Biobank ID`, `N cases`, `N controls`) %>% 
	filter((`N cases` > 100) & (`N controls` > 100))

setorder(phenotype_ancestry_information_binary , Description)

# Check the overall counts match those in the nominate candidates sheet
binary_summary <- phenotype_ancestry_information_binary %>% group_by(Description, `Phenotype ID`) %>%
	summarise(`N cases` = sum(`N cases`),
		`N controls` = sum(`N controls`), 
		`N biobanks represented` = length(unique(`Biobank ID`)),
		`N ancestries represented` = length(unique(Ancestry)))

# Binary traits
phenotype_ancestry_information_cts <- phenotype_ancestry_information %>% filter(!is.na(`N`)) %>%
	select(Description, `Phenotype ID`, Ancestry, Biobank, `Biobank ID`, `N`)
setorder(phenotype_ancestry_information_cts , Description)

# Check the overall counts match those in the nominate candidates sheet
cts_summary <- phenotype_ancestry_information_cts %>% group_by(Description, `Phenotype ID`) %>%
	summarise(`N` = sum(N),
		`N biobanks represented` = length(unique(`Biobank ID`)),
		`N ancestries represented` = length(unique(Ancestry)))

# Write to file as a csvs and import into google sheets (supplementary tables 4-7)
# Full tables
fwrite(phenotype_ancestry_information_binary, file="Tables/supplementary_table_4.tsv", sep='\t', quote=FALSE)
fwrite(phenotype_ancestry_information_cts, file="Tables/supplementary_table_5.tsv", sep='\t', quote=FALSE)

# Summary tables
fwrite(binary_summary, file="Tables/supplementary_table_6.tsv", sep='\t', quote=FALSE)
fwrite(cts_summary, file="Tables/supplementary_table_7.tsv", sep='\t', quote=FALSE)

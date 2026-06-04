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

ancestry_information <- file_info %>% group_by(dataset, ancestry) %>%
	summarise(n=max(as.integer(n), (as.integer(n_cases) + as.integer(n_controls)), na.rm=TRUE), .groups='drop')

ancestry_information <- ancestry_information %>%
	rename(
		`Biobank ID` = dataset,
		Ancestry = ancestry,
		N = n) %>%
	mutate(Biobank = unlist(renaming_plot_biobank_list[`Biobank ID`])) %>%
	select(Biobank, `Biobank ID`, Ancestry, N)

fwrite(ancestry_information, file="Tables/supplementary_table_8.tsv", sep='\t', quote=FALSE)

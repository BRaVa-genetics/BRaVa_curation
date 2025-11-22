#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("../phenotypes/BRaVa_phenotypes_utils.r")
source("meta_analysis_utils.r")

main <- function(args)
{
	data_dir <- args$analysis_results_folder
	out <- args$out

	source("meta_analysis_utils.r")
	source("../phenotypes/BRaVa_phenotypes_utils.r")

	files <- grep(".*cleaned/.*\\.gene\\..*",
		dir(data_dir, full.names=TRUE, recursive=TRUE), value=TRUE)
	dt <- rbindlist(lapply(files, extract_file_info), fill=TRUE)
	dt <- dt %>% mutate(dataset = gsub(".*/", "", dataset), file = files)

	phes <- c("AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
			"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
			"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
			"IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat",
			"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
			"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
			"AlcCons", "AST", "BMI", "CRP", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep")
	dt <- dt %>% filter(phenotype %in% phes)
	dt_list <- list()
	
	j <- 1
	if (any(duplicated(dt))) {
		cat("There are multiple matches to at least one (phenotypes, ancestry, biobank) tuple\n")
	} else {
		for (i in 1:nrow(dt)) {
			cat(paste0("determining genomic control factors for ", dt$phenotype[i],
				" in (", dt$dataset[i], ", ", dt$ancestry[i], ")\n"))
			cat(paste0("using file: ", dt$file[i], "\n"))
			print(dt[i,])
			if (!is.na(dt$n_cases[i])) {
				if ((as.integer(dt$n_cases[i]) < as.integer(args$min_cases)) | 
					(as.integer(dt$n_controls[i]) < as.integer(args$min_cases))) {
					cat("not sufficiently many cases, next trait!\n")
					next
				}
			}
							
			chisq <- qchisq(c(0.95, 0.99, 0.999), df = 1)
			dt_list[[j]] <- fread(dt$file[i]) %>% group_by(max_MAF, Group) %>% summarise(
				lambda_95_Burden = qchisq(quantile(Pvalue_Burden, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
				lambda_99_Burden = qchisq(quantile(Pvalue_Burden, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
				lambda_99.9_Burden = qchisq(quantile(Pvalue_Burden, probs=0.001, na.rm=TRUE), df=1, lower=FALSE) / chisq[3],
				lambda_95_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
				lambda_99_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
				lambda_99.9_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.001, na.rm=TRUE), df=1, lower=FALSE) / chisq[3],
				`lambda_95_SKAT-O` = qchisq(quantile(Pvalue, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
				`lambda_99_SKAT-O` = qchisq(quantile(Pvalue, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
				`lambda_99.9_SKAT-O` = qchisq(quantile(Pvalue, probs=0.001, na.rm=TRUE), df=1, lower=FALSE) / chisq[3]
			) %>% mutate(ancestry = dt$ancestry[i], dataset=dt$dataset[i],
				phenotype=dt$phenotype[i], sex=dt$sex[i])
			j <- i+1
		}
	}

	dt <- rbindlist(dt_list)
	melted_dt <- melt(dt, 
		id.vars = c("Group", "ancestry", "dataset", "phenotype", "max_MAF", "sex"), # Columns to keep
		measure.vars = c("lambda_95_Burden", "lambda_99_Burden", "lambda_99.9_Burden",
		               "lambda_95_SKAT", "lambda_99_SKAT", "lambda_99.9_SKAT",  
		               "lambda_95_SKAT-O", "lambda_99_SKAT-O", "lambda_99.9_SKAT-O"), # Columns to melt
		variable.name = "lambda_type", # Name of the new column for the melted variable names
		value.name = "lambda_value" # Name of the new column for the melted values
	)
	fwrite(melted_dt, file=args$out, sep='\t')
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks")
parser$add_argument("--out",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz",
	required=FALSE, help="Output file path")
parser$add_argument("--min_cases", required=FALSE, default=100,
	help="What is the minimum number of cases required to be included in the meta-analysis?")
args <- parser$parse_args()

main(args)
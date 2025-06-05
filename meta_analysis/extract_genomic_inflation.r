#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("../phenotypes/BRaVa_phenotypes_utils.r")
# pilot_phenotypes <- extract_BRaVa_pilot_phenotypes()

main <- function(args)
{
	data_dir <- args$analysis_results_folder
	out <- args$out
	phe <- args$phenotypeID
	
	source("meta_analysis_utils.r")
	source("../phenotypes/BRaVa_phenotypes_utils.r")

	if (is.null(phe)) {
		phes <- c(
			"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
			"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
			"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
			"IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat",
			"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
			"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
			"AlcCons", "AST", "BMI", "CRP", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep"
		)
		# phes <- pilot_phenotypes (it's the same as above)
	} else {
		phes <- phe
	}

	datasets <- c("all-of-us", "alspac", "biome", "bbj", "ckb", "ccpm",
		"decode", "egcut", "dan-rav", "genes-and-health", "gel", "pmbb",
		"mgbb", "qatar-genomes", "uk-biobank", "viking-genes")

	ancestries <- c("AFR", "AMR", "EAS", "EUR", "SAS")
	i <- 1
	dt_list <- list()
	for (phe in phes) {
		for (dataset in datasets) {
			for (anc in ancestries) {
				for (sex in c("ALL", "M", "F")) {
					file_gene <- grep(
						paste0(".*cleaned.*", dataset, "\\..*\\.",  phe, "\\..*\\.", sex, "\\..*", anc, ".*\\.gene\\..*"),
						dir(data_dir, full.names=TRUE, recursive=TRUE), value=TRUE)
					cat(paste0("determining genomic control factors for ", phe, " in (", dataset, ", ", anc, ")\n"))
					if (length(file_gene) == 1) {
						cat(paste0("using file: ", file_gene, "\n"))
						chisq <- qchisq(c(0.95, 0.99, 0.999), df = 1)
						dt_list[[i]] <- fread(file_gene) %>% group_by(max_MAF, Group) %>% summarise(
							lambda_95_Burden = qchisq(quantile(Pvalue_Burden, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
							lambda_99_Burden = qchisq(quantile(Pvalue_Burden, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
							lambda_99.9_Burden = qchisq(quantile(Pvalue_Burden, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[3],
							lambda_95_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
							lambda_99_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
							lambda_99.9_SKAT = qchisq(quantile(Pvalue_SKAT, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[3],
							lambda_95 = qchisq(quantile(Pvalue, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
							lambda_99 = qchisq(quantile(Pvalue, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
							lambda_99.9 = qchisq(quantile(Pvalue, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[3]
							) %>% mutate(ancestry = anc, dataset=dataset, phenotype=phe, sex=sex)
						i <- i+1
					}
				}
				if (length(file_gene) > 1) {
					cat("skipped: multiple matches to this (phenotypes, ancestry, biobank) tuple\n")
				}
			}
		}
	}
	dt <- rbindlist(dt_list)
	melted_dt <- melt(dt, 
		id.vars = c("Group", "ancestry", "dataset", "phenotype", "max_MAF", "sex"), # Columns to keep
		measure.vars = c("lambda_95_Burden", "lambda_99_Burden", "lambda_99.9_Burden",
		               "lambda_95_SKAT", "lambda_99_SKAT", "lambda_99.9_SKAT",  
		               "lambda_95", "lambda_99", "lambda_99.9"), # Columns to melt
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
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to evaluate. If null, this script determine lambdas for all traits")
args <- parser$parse_args()

main(args)
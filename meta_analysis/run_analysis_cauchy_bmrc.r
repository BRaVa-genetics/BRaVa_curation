#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{
	data_dir <- args$analysis_results_folder
	phe <- args$phenotypeID
	
	# Ensure that the folder is already present
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
	} else {
		phes <- phe
	}

	datasets <- c("all-of-us", "alspac", "biome", "bbj", "ckb", "ccpm",
		"decode", "egcut", "dan-rav", "genes-and-health", "gel", "pmbb",
		"mgbb", "qatar-genomes", "uk-biobank", "viking-genes")

	ancestries <- c("AFR", "AMR", "EAS", "EUR", "SAS")

	for (phe in phes) {
		for (dataset in datasets) {
			for (anc in ancestries) {
				file_gene <- grep(
					paste0(".*cleaned.*", dataset, "\\..*",  phe, ".*", anc, ".*\\.gene\\..*.txt.gz$"),
					dir(data_dir, full.names=TRUE, recursive=TRUE), value=TRUE)
				out <- gsub(".txt.gz$", ".extra_cauchy.gz", file_gene)
				cat(paste0("carrying additional Cauchy combinations for ", phe, "\n"))
				if (length(file_gene == 1)) {
					cat(paste0("using file: ", file_gene, "\n"))
					cat(paste0("output file: ", out, "\n"))
					system(paste(
						"sbatch run_cauchy_combination_of_results.sh",
						file_gene, out))
				}
				if (length(file_gene) > 1) {
					cat("skipped: multiple matches to this (phenotypes, ancestry, biobank) tuple\n")
				}
			}
		}
	}	
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to plot. If null, this script will plot everything in the folder.")
args <- parser$parse_args()

main(args)
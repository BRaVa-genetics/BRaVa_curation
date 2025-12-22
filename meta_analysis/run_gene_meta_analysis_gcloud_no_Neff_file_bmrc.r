#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("meta_analysis_utils.r")

main <- function(args)
{
	data_dir <- args$data_dir
	n_cases <- args$n_cases
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_meta_results_dir <- paste0(out_dir, "/gene/n_cases_", n_cases, "_no_Neff_file")
	
	# Ensure that the folder is already present
	system(paste("mkdir -p", out_meta_results_dir))
	source("../phenotypes/BRaVa_phenotypes_utils.r")

	# Assumes that we have the files locally within a file structure as defined in the munging scripts
	# First, let's determine the collection of phenotypes that we are testing and the 
	# collection of (population, phenotype) pairs for that biobank to include in the meta-analysis

	biobanks <- dir(data_dir)[file.info(dir(data_dir, full.names=TRUE))$isdir]
	biobanks <- intersect(biobanks, names(file_check_information$dataset))
	results_dt_list <- list()

	for (biobank in biobanks)
	{
		biobank_results_files <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"))
		biobank_results_files_full <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"), full.names=TRUE)
		results <- lapply(biobank_results_files, extract_file_info)
		results_dt_list[[biobank]] <- data.table(
			filename = biobank_results_files_full,
			phenotypeID = sapply(results, `[[`, "phenotype"),
			pop = sapply(results, `[[`, "ancestry"),
			biobank = biobank,
			sex = sapply(results, `[[`, "sex")
		)
	}

	results_dt <- rbindlist(results_dt_list)
	results_dt <- results_dt %>% filter(!(filename %in% grep("\\.extra_cauchy\\.", filename, value=TRUE)))

	if (is.null(phe)) {
		phes <- c(
			"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
			"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
			"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
			"IFHern", "ILDSarc", "MatHem", "NonRheuValv", "Pancreat",
			"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
			"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
			"AlcCons", "AST", "BMI", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep", "CRP"
		)
	} else {
		phes <- phe
	}

	# Here, we must remove any files that have been deemed to be inflated.
	dt_inflation <- fread(args$inflation_file)
	dt_inflation <- unique(dt_inflation %>% filter(Group == "synonymous") %>% 
		filter(max_MAF != 0.01,
			lambda_value > 1.3,
			!(lambda_type %in% c("lambda_50_Burden", "lambda_50_SKAT", "lambda_50"))) %>% 
		select(phenotype, dataset, ancestry, sex))
	# Manual curation, adding the following (biobank, trait) tuples containing spurious 
	# associations
	dt_inflation <- rbind(dt_inflation, data.table(
		phenotype = c("ColonRectCanc", "Height"),
		dataset = c("egcut", "mgbb"),
		ancestry = c("EUR", "AMR"),
		sex = c("ALL", "ALL"))
	) %>% rename(phenotypeID = phenotype, biobank = dataset, pop = ancestry)
	dt_inflation <- setdiff(dt_inflation, data.table(
	phenotypeID = c("Height"),
	biobank = c("uk-biobank"),
	pop = c("EUR"),
	sex = c("ALL")))

	# Remove those (phenotype, biobank, sex) files from the meta-analysis
	setkeyv(dt_inflation, c("phenotypeID", "biobank", "pop", "sex"))
	setkeyv(results_dt,  c("phenotypeID", "biobank", "pop", "sex"))
	results_dt <- setdiff(results_dt, merge(dt_inflation, results_dt))

	# Everything
	for (phe in phes) {
		for (s in c("ALL", "M", "F")) {
			files_gene <- (results_dt %>% filter(phenotypeID == phe, sex == s))$filename
			if (length(files_gene) <= 1) { 
				cat("Either the phenotype is not present, or there is only a single file for:\n")
				cat(phe, s, "\n")
			} else {
				files_gene <- paste(files_gene, collapse=",")
				out <- paste0(out_meta_results_dir, "/", phe, "_", s, "_gene_meta_analysis_", n_cases, "_cutoff_no_Neff_file.tsv.gz")
				cat(paste0("carrying out meta-analysis of ", phe, " in ", s, "\n"))
				cat(paste0("\nFiles in the analysis: ",
					paste0(strsplit(files_gene, split=",")[[1]], collapse='\n'), "\n"))
				system(paste(
					"sbatch run_meta_analysis_gcloud_no_Neff_file_bmrc.sh",
					files_gene, out))
				cat(paste0("submitted meta-analysis of ", phe, " completed\n\n"))
			}
		}
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--data_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks",
	required=FALSE, help="Location of the folder containing the input into the meta-analysis")
parser$add_argument("--n_cases", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--out_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs",
	required=FALSE, help="Output folder path")
parser$add_argument("--inflation_file", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz",
	required=FALSE, help="Inflation file")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to run meta-analysis on. Note: must match the naming in input folder.")
args <- parser$parse_args()

main(args)

#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{
	data_dir <- args$analysis_results_folder
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_plot_dir <- paste0(out_dir, "/gene")
	
	# Ensure that the folder is already present
	system(paste("mkdir -p", out_plot_dir))
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
			"AlcCons", "AST", "BMI", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep", "CRP"
		)
	} else {
		phes <- phe
	}

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
			sex = sapply(results, `[[`, "sex"),
			last_name = sapply(results, `[[`, "last_name")
		)
	}

	results_dt <- rbindlist(results_dt_list)
	results_dt <- results_dt %>% filter(!(filename %in% grep("\\.extra_cauchy\\.", filename, value=TRUE)))
	results_dt <- results_dt %>% filter(phenotypeID %in% phes)

	# Ensure that there are not multiple results files for a given (phenotype, biobank, sex) combination
	if (nrow(results_dt %>% group_by(phenotypeID, pop, biobank, sex, last_name) %>% filter(n() > 1)) > 0) {
		message("Multiple (pheno, pop, biobank, sex, last_name) files - this will be a problem in meta-analysis")
		quit(status = 1)
	}

	# Everything
	for (i in 1:nrow(results_dt)) {
		row <- results_dt[i,]
		filename <- row$filename
		out <- paste0(out_plot_dir, "/", row$biobank, "_", row$phenotypeID, "_",
			row$sex, "_", row$pop, "_", row$last_name, "_gene_qq.pdf")
		cat(paste0("using file: ", filename, "\n"))
		cmd <- paste("sbatch run_analysis_qq_gcloud_bmrc.sh", filename, out)
		system(cmd)
		cat(paste0(cmd, "\n"))
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks")
parser$add_argument("--out_dir",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plots_biobank_specific",
	required=FALSE, help="Output folder path")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to plot. If null, this script will plot everything in the folder.")
args <- parser$parse_args()

main(args)

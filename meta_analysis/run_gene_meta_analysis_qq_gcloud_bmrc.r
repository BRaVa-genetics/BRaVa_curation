#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{
	data_dir <- args$meta_analysis_results_folder
	n_cases <- args$n_cases
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_plot_dir <- paste0(out_dir, "/gene/n_cases_", n_cases)
	
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
			"AlcCons", "AST", "BMI", "CRP", "HDLC", "Height", "LDLC",
			"TChol", "TG", "WHRBMI", "HipRep"
		)
	} else {
		phes <- phe
	}

	for (phe in phes) {
		file_gene <- grep(phe, dir(data_dir, full.names=TRUE), value=TRUE)
		if (any(grepl("extra_cauchy", file_gene))) {
			file_gene <- file_gene[-grep("extra_cauchy", file_gene)]
		}
		for (s in c("ALL", "M", "F")) {
			if (sum(grepl(paste0("_", s, "_"), file_gene)) == 1) {
				out <- paste0(out_plot_dir, "/", phe, "_", s, "_gene_meta_analysis_qq.pdf")
				file_gene_tmp <- file_gene[which(grepl(paste0("_", s, "_"), file_gene))]
				cat(paste0("carrying out plotting of gene QQ for ", phe, " in ", s, "\n"))
				cat(paste0("using file: ", file_gene_tmp, "\n"))
				system(paste(
					"sbatch run_meta_analysis_qq_gcloud_bmrc.sh",
					file_gene_tmp, out))
			}
		}
		cat(paste0("submitted meta-analysis QQ plotting of ", phe, "\n\n"))
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--meta_analysis_results_folder",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	required=FALSE)
parser$add_argument("--out_dir",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plots",
	required=FALSE, help="Output folder path")
parser$add_argument("--n_cases", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to plot. If null, this script will plot everything in the folder.")
args <- parser$parse_args()

main(args)

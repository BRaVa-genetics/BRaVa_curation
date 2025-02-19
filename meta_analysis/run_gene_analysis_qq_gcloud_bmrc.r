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
		cat(paste0("carrying out plotting of gene QQ for ", phe, "\n"))
		for (dataset in datasets) {
			for (anc in ancestries) {
				file_gene <- grep(
					paste0(".*cleaned.*", dataset, "\\..*",  phe, ".*", anc),
					dir(data_dir, full.names=TRUE, recursive=TRUE), value=TRUE)
				if ((length(file_gene) > 0) & (length(file_gene) <= 3)) {
					sexes <- c()
					for (f in file_gene) {
						sexes <- c(sexes, ifelse(
							grepl("\\.ALL\\.", f), "ALL",
								ifelse(grepl("\\.F\\.", f), "F",
									ifelse(grepl("\\.M\\.", f), "M", NA))))
					}
					if (length(unique(sexes)) == length(sexes)) {
						# output file
						out <- paste0(out_plot_dir, "/", dataset, "_", phe, "_",
							sexes, "_", anc, "_gene_meta_analysis_qq.pdf")
						for (i in 1:length(file_gene)) {
							cat(paste0("using file: ", file_gene[i], "\n"))
							cmd <- paste("sbatch run_analysis_qq_gcloud_bmrc.sh",
								file_genep[i], out[i])
							system(cmd)
							cat(paste0(cmd, "\n"))
							cat(paste0("submitted meta-analysis QQ plotting of ", phe, "\n\n"))
						}
					} else {
						cat(paste0("There are multiple files for one of the sexes for (",
							phe, ", ", dataset, ", ", anc, ")\n"))
						print(file_gene)
					}
				} else {
					if (length(file_gene) == 0) {
						cat(paste0("This phenotype is not available for (",
							phe, ", ", dataset, ", ", anc, ")\n"))
					} else {
						cat(paste0("There are more than 3 unique files for (",
							phe, ", ", dataset, ", ", anc, ")\n"))
					}
				}
			}
		}
	}	
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs")
parser$add_argument("--out_dir",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plots_biobank_specific",
	required=FALSE, help="Output folder path")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to plot. If null, this script will plot everything in the folder.")
args <- parser$parse_args()

main(args)

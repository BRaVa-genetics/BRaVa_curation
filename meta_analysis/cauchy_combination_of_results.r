#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

source("meta_analysis_utils.r")

main <- function(args)
{
	# Create a list to loop over
	# Throw an error if the number of entries of the passed vectors are different
	max_MAF <- strsplit(args$max_MAF_list, split=":")[[1]]
	cauchy_n <- length(max_MAF)

	group <- strsplit(args$Group_list, split=":")[[1]]
	names <- strsplit(args$names_list, split=":")[[1]]

	if (!all(c(length(group), length(names)) == cauchy_n)) {
		stop("the list lengths passed do not match each other")
	}

	cauchy_list <- list()
	for (i in 1:cauchy_n) {
		cauchy_list[[names[i]]] <- list(
			max_MAF = strsplit(max_MAF[i], split=",")[[1]],
			Group = strsplit(group[i], split=",")[[1]]
		)
	}

	# Now, read in the results file and apply the Cauchy combination to it
	dt <- fread(args$results_file)
	dt_cauchy_list <- list()

	for (n in names(cauchy_list)) {
		dt_cauchy_list[[n]] <- dt %>% filter(
			Group %in% cauchy_list[[n]][["Group"]],
			max_MAF %in% cauchy_list[[n]][["max_MAF"]],
			!is.na(Pvalue)) %>% 
		group_by(Region) %>% 
		summarise(
			Pvalue = cauchy_combination(Pvalue),
			Pvalue_Burden = cauchy_combination(Pvalue_Burden),
			Pvalue_SKAT = cauchy_combination(Pvalue_SKAT),
			number_of_pvals = n()
		) %>% mutate(
		Pvalue = ifelse(Pvalue > 1e+15,
			(1/Pvalue) / pi,
			pcauchy(Pvalue, lower.tail=FALSE)),
		Pvalue_Burden = ifelse(Pvalue_Burden > 1e+15,
			(1/Pvalue_Burden) / pi,
			pcauchy(Pvalue_Burden, lower.tail=FALSE)),
		Pvalue_SKAT = ifelse(Pvalue_SKAT > 1e+15,
			(1/Pvalue_SKAT) / pi,
			pcauchy(Pvalue_SKAT, lower.tail=FALSE))) %>% 
		mutate(
			Pvalue = ifelse(Pvalue > (1 - 1e-10),
				(1 - 1/number_of_pvals), Pvalue),
			Pvalue_Burden = ifelse(Pvalue_Burden > (1 - 1e-10),
				(1 - 1/number_of_pvals), Pvalue_Burden),
			Pvalue_SKAT = ifelse(Pvalue_SKAT > (1 - 1e-10),
				(1 - 1/number_of_pvals), Pvalue_SKAT)) %>% 
		select(-number_of_pvals) %>% mutate(Group = n)
	}
	dt_cauchy <- rbindlist(dt_cauchy_list)
	fwrite(dt_cauchy, file=args$out, sep='\t')
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--results_file",
    default=NULL, required=FALSE,
    help="File containing the summary statistics for Cauchy combination")
parser$add_argument("--max_MAF_list",
    default="0.001,1e-04:0.01,0.001,1e-04", required=FALSE,
    help="List of vectors of MAFs for inclusion in the Cauchy combination")
parser$add_argument("--Group_list",
    default="pLoF,damaging_missense_or_protein_altering,pLoF;damaging_missense_or_protein_altering:pLoF,damaging_missense_or_protein_altering,pLoF;damaging_missense_or_protein_altering", required=FALSE,
    help="List of vectors of Groups for inclusion in the Cauchy combination")
parser$add_argument("--names_list",
    default="pLoF_damaging_0.001:pLoF_damaging_0.01", required=FALSE,
    help="Meta-analysis results file")
parser$add_argument("--out",
    default="cauchy_combination.tsv.gz",
    required=FALSE, help="Output filepath")
args <- parser$parse_args()

main(args)

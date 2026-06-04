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
	files <- grep(".tsv.gz$", dir(data_dir, full.names=TRUE), value=TRUE)
	chisq <- qchisq(c(0.95, 0.99, 0.999), df = 1)

	dt_list <- list()
	for (file_gene in files) {
		cat(paste0("determining genomic control factors for ", basename(file_gene)))
		dt_list[[file_gene]] <- fread(file_gene) %>% group_by(max_MAF, Group, type, class) %>% 
		summarise(
			lambda_95 = qchisq(quantile(Pvalue, probs=0.05, na.rm=TRUE), df=1, lower=FALSE) / chisq[1],
			lambda_99 = qchisq(quantile(Pvalue, probs=0.01, na.rm=TRUE), df=1, lower=FALSE) / chisq[2],
			lambda_99.9 = qchisq(quantile(Pvalue, probs=0.001, na.rm=TRUE), df=1, lower=FALSE) / chisq[3],
		) %>% mutate(phenotype=gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file_gene))
	}
	dt <- rbindlist(dt_list)
	melted_dt <- melt(dt, 
		id.vars = c("Group", "phenotype", "max_MAF", "class", "type"), # Columns to keep
		measure.vars = c("lambda_95", "lambda_99", "lambda_99.9"),
		variable.name = "lambda_type", # Name of the new column for the melted variable names
		value.name = "lambda_value" # Name of the new column for the melted values
	) %>% mutate(lambda_type = paste(lambda_type, class, sep="_"))
	fwrite(melted_dt, file=args$out, sep='\t')
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder", required=FALSE,
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100")
parser$add_argument("--out",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries_meta.tsv.gz",
	required=FALSE, help="Output file path")
args <- parser$parse_args()

main(args)
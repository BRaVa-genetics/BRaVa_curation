#!/bin/Rscript
library(data.table)
library(dplyr)

# Names of the phenotypes to run
data_dir <- "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud"
n_cases <- 100
out_meta_results_dir <- paste0(
	"~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/n_cases_", n_cases, "_nglmm")
system(paste("mkdir -p", out_meta_results_dir))
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")
source("~/Repositories/BRaVa_curation/phenotypes/BRaVa_phenotypes_utils.r")

# Assumes that we have the files locally within a file structure as defined in the munging scripts
# First, let's determine the collection of phenotypes that we are testing and the 
# collection of (population, phenotype) pairs for that biobank to include in the meta-analysis

biobanks <- dir(data_dir)
results_dt_list <- list()

for (biobank in biobanks)
{
	biobank_results_files <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"))
	biobank_results_files_full <- dir(paste0(data_dir, "/", biobank, "/cleaned/gene"), full.names=TRUE)
	results <- lapply(biobank_results_files, extract_file_info)
	results_dt_list[[biobank]] <- data.table(filename = biobank_results_files_full, phenotypeID = sapply(results, `[[`, 4), pop = sapply(results, `[[`, 7))
}

results_dt <- rbindlist(results_dt_list)
if (length(grep("\\.extra_cauchy\\.", results_dt$filename)) > 0) {
	results_dt <- results_dt %>% filter(-grep("\\.extra_cauchy\\.", filename))
}
# pilot_phenotypes <- extract_BRaVa_pilot_phenotypes()
# phenotypeIDs <- intersect(pilot_phenotypes, unique(results_dt$phenotypeID))
phenotypeIDs <- "AMD"

for (phe in phenotypeIDs)
{
	files_gene <- (results_dt %>% filter(phenotypeID == phe))$filename
	files_gene <- paste(files_gene, collapse=",")
	cat(paste0("carrying out meta-analysis of ", phe, "\n"))
	cat(paste0("\nFiles in the analysis: ", paste0(strsplit(files_gene, split=",")[[1]], collapse='\n')))
	system(paste(
		"Rscript meta_analysis.r",
		"--file_paths", files_gene,
		"--no_sex_check ",
		"--Neff_weights_file", "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/Neff_weights.tsv.gz",
		"--out", paste0(
			out_meta_results_dir, "/", phe, "_gene_meta_analysis_", n_cases, "_cutoff.tsv.gz")
		)
	)
	cat(paste0("meta-analysis of ", phe, " completed\n\n"))
	break
}

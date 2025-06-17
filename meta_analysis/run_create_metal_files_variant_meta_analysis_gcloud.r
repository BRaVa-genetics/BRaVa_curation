#!/bin/Rscript
library(data.table)
library(dplyr)

# Names of the phenotypes to run
data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks")
out_meta_results_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant")

system(paste("mkdir -p", out_meta_results_dir))
source("meta_analysis_utils.r")
source("../phenotypes/BRaVa_phenotypes_utils.r")

# Assumes that we have the files locally within a file structure as defined in the munging scripts
# First, let's determine the collection of phenotypes that we are testing and the 
# collection of (population, phenotype) pairs for that biobank to include in the meta-analysis

biobanks <- dir(data_dir)
results_dt_list <- list()
min_cases <- 100

for (biobank in biobanks)
{
	biobank_results_files <- dir(paste0(data_dir, "/", biobank, "/cleaned/variant"))
	biobank_results_files_full <- dir(paste0(data_dir, "/", biobank, "/cleaned/variant"), full.names=TRUE)
	results <- lapply(biobank_results_files, extract_file_info)
	results_dt_list[[biobank]] <- data.table(
		filename = biobank_results_files_full,
		phenotypeID = sapply(results, `[`, 4),
		pop = sapply(results, `[`, 7),
		binary = sapply(results, `[`, "binary")
		) %>% 
	mutate(n_cases = ifelse(binary, sapply(results, `[`, "n_cases"), NA),
		n_controls = ifelse(binary, sapply(results, `[`, "n_controls"), NA),
		n = ifelse(binary, NA, sapply(results, `[`, "n"))) %>% 
	filter((n > min_cases) | ((n_controls > min_cases) & (n_cases > min_cases)))
}

results_dt <- rbindlist(results_dt_list)
phenotypeIDs <- intersect(BRaVa_pilot_phenotypes, unique(unlist(results_dt$phenotypeID)))

for (phe in phenotypeIDs)
{
	files_variant <- (results_dt %>% filter(phenotypeID == phe))$filename
	files_variant <- paste(files_variant, collapse=",")

	cat(paste0("carrying out meta-analysis of ", phe, "\n"))
	cat(paste0("creating METAL file...\n"))

	system(paste(
		"Rscript create_metal_script.r ",
		"--phenotypeID", phe,
		"--files", files_variant,
		"--out_folder", out_meta_results_dir
		)
	)
	cat(paste0("METAL file creation for meta-analysis of ", phe, " completed\n\n"))
}

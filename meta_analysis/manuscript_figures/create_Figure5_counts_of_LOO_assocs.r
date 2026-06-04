#!/bin/Rscript
library(data.table)
library(dplyr)

source("../meta_analysis_utils.r")

files_list <- c(
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_all-of-us",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_bbj",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_biome",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_ccpm", 
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_egcut",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_gel", 
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_genes-and-health",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_mgbb",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_pmbb",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/minus_uk-biobank")

meta_list <- list()
for (file in files_list)
{
	files <- dir(file, full.names=TRUE)
	files <- grep("tsv.gz", files, value=TRUE)
	meta_list[[file]] <- list()
	for (f in files) {
		cat(f, "\n")
		phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", f)
		meta_list[[file]][[f]] <- fread(f) %>% filter(
			max_MAF %in% c("1e-04", "0.001"),
			Group %in% c(
				"damaging_missense_or_protein_altering",
				"pLoF",
				"pLoF;damaging_missense_or_protein_altering"),
			((class == "Burden" & type == "Inverse variance weighted") |
			 (class %in% c("SKAT", "SKAT-O") & type == "Stouffer"))) %>%
		mutate(phenotype = phenotype) %>% mutate(hit = Pvalue < ((0.05/20000)/(2*3*3)))
	}
	meta_list[[file]] <- rbindlist(meta_list[[file]])
}

for (biobank in names(meta_list)) {
	meta_list[[biobank]]$biobank <- ifelse(biobank == "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
		"all", gsub(".*/minus_", "", biobank))
}
meta <- rbindlist(meta_list)

fwrite(meta %>% filter(hit), sep = "\t", quote=FALSE,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_all_meta_subsets_101025_LOO.tsv.gz")

meta <- meta %>% filter(!is.na(Pvalue)) %>%
	mutate(Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue), weights = 1)
meta_cauchy <- run_cauchy(meta %>% group_by(Region, phenotype, biobank), "weights", "Cauchy_stat", "Pvalue", "Cauchy_Pvalue")
meta_cauchy <- meta_cauchy %>% mutate(hit = Cauchy_Pvalue < 2.5e-6)

fwrite(meta_cauchy %>% filter(hit), sep = "\t", quote=FALSE,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_cauchy_assocs_from_all_meta_subsets_101025_LOO.tsv.gz")

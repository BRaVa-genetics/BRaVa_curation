#!/bin/Rscript
library(data.table)
library(dplyr)

source("../meta_analysis_utils.r")

files_list <- c("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/just_uk-biobank_and_all-of-us",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AMR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AFR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EUR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/SAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/non_EUR")

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
				"pLoF;damaging_missense_or_protein_altering")) %>%
		mutate(phenotype = phenotype) %>% filter(
			(Pvalue < 6.7e-7 & class == "Burden" & type == "Inverse variance weighted") |
			(Pvalue < 2.5e-7 & class %in% c("SKAT", "SKAT-O") & type == "Stouffer"))
	}
	meta_list[[file]] <- rbindlist(meta_list[[file]])
}

for (anc in names(meta_list)) {
	meta_list[[anc]]$ancestry <- ifelse(anc == "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
		"all", gsub(".*/", "", anc))
}
meta <- rbindlist(meta_list)

fwrite(meta, sep = "\t", quote=FALSE,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_all_meta_subsets_101025.tsv.gz")

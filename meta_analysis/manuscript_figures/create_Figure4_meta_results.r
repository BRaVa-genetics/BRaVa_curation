#!/bin/Rscript

# Create the manhattan plots
# What do we need for this?
# The region, the gene name, the minimum P-value among the hits

# Install and load biomaRt
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(ggplot2)
devtools::install_github("mkanai/rgsutil")
library(rgsutil)
source("../meta_analysis_utils.r")
source("../QC/utils/pretty_plotting.r")

# Meta-analysis files for plotting should have been copied from BMRC to gcloud
# gsutil cp -r gene gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/
# gsutil mv 'gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/*gz' gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/meta-analysis/

# Connect to Ensembl BioMart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
cloud <- FALSE

# Choose the dataset for human genes
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Retrieve Ensembl IDs, start and end positions
gene_info <- getBM(attributes = c(
	"ensembl_gene_id",
	"external_gene_name",
	"chromosome_name",
	"start_position",
	"end_position"),
	mart = ensembl_dataset)

# Display the first few rows of the result
gene_info <- data.table(gene_info, key = "ensembl_gene_id")

# Download the results files using Masa's rgsutil
# file_paths <- c(
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/meta-analysis/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/AFR/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/AMR/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/EAS/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/EUR/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/SAS/",
# 	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/non_EUR/"
# )

file_paths <- c(
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AFR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AMR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EUR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/SAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/non_EUR"
)

file_root <- c("meta_analysis", "AFR", "AMR", "EAS", "EUR", "SAS", "non_EUR")

for (i in 1:length(file_paths)) {
	# First, grab the files
	if (cloud) {
		files <- system(paste("gsutil ls", file_paths[i]), intern=TRUE)
	} else {
		files <- grep(".gz$", dir(file_paths[i], full.names=TRUE), value=TRUE)
	}
	meta_list <- list()

	for (file in files) {
		cat(file, "\n")
		phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", file)
		if (cloud) {
			meta_list[[file]] <- rgsutil::read_gsfile(file)
		} else {
			meta_list[[file]] <- fread(file)
		}
		meta_list[[file]] <- meta_list[[file]] %>% filter(
			max_MAF %in% c("1e-04", "0.001"),
			Group %in% c(
				"damaging_missense_or_protein_altering",
				"pLoF",
				"pLoF;damaging_missense_or_protein_altering")) %>%
		filter(
			((class == "Burden") & (type == "Inverse variance weighted")) | 
			((class != "Burden") & (type == "Stouffer"))) %>%
		mutate(phenotype = phenotype) %>% group_by(Region, phenotype) %>% filter(Pvalue == min(Pvalue))
	}

	meta_list <- rbindlist(meta_list) %>% rename(ensembl_gene_id = Region)
	setkey(meta_list, "ensembl_gene_id")
	meta_list <- merge(gene_info, meta_list)
	meta_list <- meta_list %>% filter(chromosome_name %in% c(seq(1,22), "X"))

	# Remove the mosaic genes:
	meta_list <- meta_list %>% filter(!ensembl_gene_id %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))
	meta_list <- meta_list %>% mutate(Pvalue = ifelse(Pvalue == 0, 1e-320, Pvalue))
	meta_list <- meta_list %>% mutate(phenotype = gsub("_.*", "", phenotype))
	meta_list <- meta_list %>% mutate(
		case_control = ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

	# Write the information to disk, to speed up recreation of the plots.
	fwrite(meta_list, quote=FALSE, sep='\t',
		file=paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/", file_root[i], "_figure_4.tsv.gz"))
}

#!/bin/Rscript
library(data.table)
library(dplyr)
library(biomaRt)

# Plan is to loop through all the results files and provide a superset of candidates to run conditional analysis with
# It's easier to run a superset than ask for more later on...

files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100", full.names=TRUE, recursive=TRUE)
files <- files[-grep("/minus_", files)]

extract_file_info_meta <- function(filename)
{
	gz <- ifelse(grepl(".gz$", filename), TRUE, FALSE)
	filename <- gsub(".gz$", "", filename)
	filename <- gsub("cleaned.", "", filename)
	file_info <- as.list(strsplit(filename, split="\\.")[[1]])
	file_info[[1]] <- gsub(".*/", "", file_info[[1]])
	file_info$phenotype <- gsub("^([A-Za-z0-9]+).*", "\\1", file_info[[1]])
	file_info$meta <- ifelse(file_info[[2]] == "tsv", "ALL", file_info[[2]])
	file_info$gz <- gz
	file_info <- file_info[c("phenotype", "meta", "gz")]
	return(file_info)
}

P <- 6.7e-7
dt_list <- list()
for (file in files) {
	cat(paste0(file, "\n"))
	# Read in the file, determine if there are any experiment-wise significant associations
	# Note that we only care about associations due to either, pLoF, pLoF:damaging_missense_or_protein_altering,
	# damaging_missense_or_protein_altering
	# We also only consider inverse variance weighting for the burden tests, and Stouffer for SKAT and SKAT-O.
	dt <- unique(fread(file) %>% 
	filter(
		Pvalue < P,
		Group %in% c(
			"pLoF",
			"damaging_missense_or_protein_altering",
			"pLoF;damaging_missense_or_protein_altering")) %>% 
	mutate(max_MAF = as.numeric(max_MAF)) %>%
	filter(max_MAF %in% c(0.001, 0.0001)) %>% 
	filter(((type == "Stouffer") & (class %in% c("SKAT", "SKAT-O"))) |
		(type == "Inverse variance weighted") & (class == "Burden")) %>%
	dplyr::select(Region))
	file_info <- extract_file_info_meta(file)
	dt$phenotype <- file_info$phenotype
	dt$meta <- file_info$meta
	dt$file <- file
	dt_list[[file]] <- dt
}

dt_final <- rbindlist(dt_list)
fwrite(dt_final, "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/candidates_from_meta_101025.tsv.gz", sep='\t', quote=FALSE)

# Double check to ensure that the ALL meta results are all present within this less stringent superset
dt_all <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_full_meta_060625.tsv.gz")
dt_all_unique <- unique(dt_all %>% dplyr::select(Region, phenotype) %>% mutate(phenotype = gsub("_.*", "", phenotype)))

# Then filter to the subset of unique phenotype gene rows
dt_final <- unique(dt_final %>% dplyr::select(phenotype, Region))

setdiff(dt_final, dt_all_unique)
setdiff(dt_all_unique, dt_final) # This should be empty

# Connect to Ensembl BioMart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host="https://asia.ensembl.org")

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
gene_info_list <- list()
gene_info_list[gene_info$ensembl_gene_id] <- gene_info$external_gene_name
dt_final$external_gene_name <- gene_info_list[dt_final$Region]
fwrite(dt_final,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_101025.csv.gz",
	sep=",", quote=FALSE)

#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(biomaRt)

source("../meta_analysis_utils.r")

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

# Count the number of hits, split by all of the data sets and ancestries
files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks",
	full.names=TRUE, recursive=TRUE)
files <- grep("/cleaned/", files, value=TRUE)
gene_files <- grep("\\.gene\\.", files, value=TRUE)

file_info <- rbindlist(lapply(gene_files, extract_file_info), fill=TRUE)
file_info <- file_info %>% mutate(filename = gene_files,
	dataset = gsub(".*/", "", dataset))

# Munge the data into a format that can be used by run_cauchy.
groups_for_cauchy <- c("pLoF", "damaging_missense_or_protein_altering",
	"pLoF;damaging_missense_or_protein_altering")
max_MAFs_for_cauchy <- c(1e-4, 1e-3)

extract_cauchy <- function(filename, Cauchy_cutoff = 2.5e-6) {
	dt <- fread(filename) %>% filter(Group %in% groups_for_cauchy, max_MAF %in% max_MAFs_for_cauchy)
	melted_dt <- rbind(
		dt %>% dplyr::select(-c("Pvalue_SKAT", "Pvalue_Burden")) %>% mutate(Test = "Pvalue"),
		dt %>% dplyr::select(-c("Pvalue", "Pvalue_Burden")) %>% rename(Pvalue = Pvalue_SKAT) %>% 
			mutate(Test = "Pvalue_SKAT"),
		dt %>% dplyr::select(-c("Pvalue_SKAT", "Pvalue")) %>% rename(Pvalue = Pvalue_Burden) %>% 
			mutate(Test = "Pvalue_Burden")
		) %>% mutate(weights = 1)
	melted_dt <- melted_dt %>% filter(!is.na(Pvalue)) %>% 
		mutate(Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue)) %>% group_by(Region)
	dt <- run_cauchy(melted_dt, "weights", "Cauchy_stat", "Pvalue", "Cauchy_Pvalue") %>%
		filter(Cauchy_Pvalue < Cauchy_cutoff)
	return(dt)
}

extract_hits <- function(filename, P_SKAT_cutoff=(0.05/20000)/(3*3*2),
	P_SKAT_O_cutoff=(0.05/20000)/(3*3*2), P_Burden_cutoff=(0.05/20000)/(3*3*2)) {
	dt <- fread(filename) %>% filter(
		(Pvalue < P_SKAT_O_cutoff) | 
		(Pvalue_SKAT < P_SKAT_cutoff) |
		(Pvalue_Burden < P_Burden_cutoff))
	melted_dt <- rbind(
		dt %>% dplyr::select(-c("Pvalue_SKAT", "Pvalue_Burden")) %>% 
			mutate(Test = "Pvalue") %>% filter(Pvalue < P_SKAT_O_cutoff),
		dt %>% dplyr::select(-c("Pvalue", "Pvalue_Burden")) %>%
			rename(Pvalue = Pvalue_SKAT) %>% mutate(Test = "Pvalue_SKAT") %>%
			filter(Pvalue < P_SKAT_cutoff),
		dt %>% dplyr::select(-c("Pvalue_SKAT", "Pvalue")) %>%
			rename(Pvalue = Pvalue_Burden) %>% mutate(Test = "Pvalue_Burden") %>%
			filter(Pvalue < P_Burden_cutoff)
		)
	return(melted_dt)
}

# Bonferroni
dt_gene_hits_all_list <- list()
for (i in 1:nrow(file_info)) {
	cat(paste0(file_info$filename[i], "\n",
		file_info$dataset[i], "\n",
		file_info$phenotype[i], "\n\n"))
	melted_dt <- extract_hits(file_info$filename[i])
	melted_dt <- melted_dt %>% 
		mutate(
			dataset = file_info$dataset[i],
			phenotype = file_info$phenotype[i],
			ancestry = file_info$ancestry[i],
			sex = file_info$sex[i]
		)
	dt_gene_hits_all_list[[i]] <- melted_dt
}

dt_gene_hits_all <- rbindlist(dt_gene_hits_all_list, fill=TRUE)
dt_gene_hits_all <- dt_gene_hits_all %>% filter(max_MAF %in% c("0.001", "1e-04"))
setkeyv(dt_gene_hits_all, c("phenotype", "dataset", "ancestry", "sex"))

# Cauchy
dt_gene_cauchy_hits_all_list <- list()
for (i in 1:nrow(file_info)) {
	cat(paste0(file_info$filename[i], "\n",
		file_info$dataset[i], "\n",
		file_info$phenotype[i], "\n\n"))
	melted_dt <- extract_cauchy(file_info$filename[i])
	melted_dt <- melted_dt %>%
		mutate(
			dataset = file_info$dataset[i],
			phenotype = file_info$phenotype[i],
			ancestry = file_info$ancestry[i],
			sex = file_info$sex[i]
		)
	dt_gene_cauchy_hits_all_list[[i]] <- melted_dt
}

dt_gene_cauchy_hits_all <- rbindlist(dt_gene_cauchy_hits_all_list, fill=TRUE)
setkeyv(dt_gene_cauchy_hits_all, c("phenotype", "dataset", "ancestry", "sex"))

# Code to generate inflation_summaries.tsv.gz is 'extract_genomic_inflation.r' in
# the folder above.
# Run the code on everything using Rscript extract_genomic_inflation.r to update it.
# Here, we must remove any files that have been deemed to be inflated.
dt_inflation <- fread(
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz")
dt_inflation <- unique(dt_inflation %>% filter(Group == "synonymous") %>% 
	filter(max_MAF != 0.01, lambda_value > 1.3) %>% 
	select(phenotype, dataset, ancestry, sex))
# Manual curation, adding the following (biobank, trait) tuples containing spurious 
# associations
dt_inflation <- rbind(dt_inflation, data.table(
	phenotype = c("ColonRectCanc", "Height"),
	dataset = c("egcut", "mgbb"),
	ancestry = c("EUR", "AMR"),
	sex = c("ALL", "ALL")))
dt_inflation <- setdiff(dt_inflation, data.table(
	phenotype = c("Height"),
	dataset = c("uk-biobank"),
	ancestry = c("EUR"),
	sex = c("ALL")))

# Remove those (phenotype, biobank, sex) files from the meta-analysis
setkeyv(dt_inflation, c("phenotype", "dataset", "ancestry", "sex"))
dt_gene_hits_all <- setdiff(dt_gene_hits_all, merge(dt_gene_hits_all, dt_inflation))
dt_gene_cauchy_hits_all <- setdiff(dt_gene_cauchy_hits_all, merge(dt_gene_cauchy_hits_all, dt_inflation))

dt_gene_hits_all <- dt_gene_hits_all %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

dt_gene_cauchy_hits_all <- dt_gene_cauchy_hits_all %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

# Write these results
fwrite(dt_gene_hits_all,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_in_all_biobanks.tsv.gz",
	sep='\t', quote=FALSE)
fwrite(dt_gene_cauchy_hits_all,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_cauchy_assocs_in_all_biobanks.tsv.gz",
	sep='\t', quote=FALSE) # This is never used (yet) - but useful to have!

# This is the meta-analysis results
files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	full.names=TRUE)
files <- grep("cutoff.tsv.gz", files, value=TRUE)

meta_list <- list()
for (file in files) {
	cat(file, "\n")
	phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", file)
	meta_list[[file]] <- fread(file) %>% filter(
		max_MAF %in% c("1e-04", "0.001"),
		Group %in% c(
			"damaging_missense_or_protein_altering",
			"pLoF",
			"pLoF;damaging_missense_or_protein_altering")) %>%
	mutate(phenotype = phenotype) %>% mutate(hit = (
		(Pvalue < ((0.05/20000)/(2*3*3)) & class == "Burden" & type == "Inverse variance weighted") |
		(Pvalue < ((0.05/20000)/(2*3*3)) & class %in% c("SKAT", "SKAT-O") & type == "Stouffer")))
}

meta_list <- rbindlist(meta_list) %>% mutate(case_control =
	ifelse(gsub("_.*", "", phenotype) %in% case_ctrl, TRUE,
		ifelse(gsub("_.*", "", phenotype) %in% cts, FALSE, NA)))

meta_list <- meta_list %>% filter(!is.na(Pvalue)) %>% 
	filter((class == "Burden" & type == "Inverse variance weighted") |
		(class %in% c("SKAT", "SKAT-O") & type == "Stouffer")) %>%
	mutate(Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue), weights = 1)
meta_cauchy <- run_cauchy(meta_list %>% group_by(Region, phenotype), "weights", "Cauchy_stat", "Pvalue", "Cauchy_Pvalue")
meta_cauchy <- meta_cauchy %>% mutate(hit = (Cauchy_Pvalue < 2.5e-6))

meta_list <- meta_list %>% rename(ensembl_gene_id = Region)
setkey(meta_list, "ensembl_gene_id")
meta_list <- merge(gene_info, meta_list, all.y=TRUE)

meta_cauchy <- meta_cauchy %>% rename(ensembl_gene_id = Region)
setDT(meta_cauchy)
setkey(meta_cauchy, "ensembl_gene_id")
meta_cauchy <- merge(gene_info, meta_cauchy, all.y=TRUE)

# Write the results
fwrite(file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_full_meta_051125.tsv.gz",
	meta_list %>% filter(hit), sep='\t', quote=FALSE)
fwrite(file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_cauchy_assocs_from_full_meta_051125.tsv.gz",
	meta_cauchy %>% filter(hit), sep='\t', quote=FALSE)

extract_hit_counts_for_figure <- function(biobank_hits, meta) {

	meta_unique <- meta %>% filter(hit) %>% group_by(case_control) %>% 
		filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
		summarise(count = length(unique(paste(Region, phenotype))))
	meta_unique_no_height <- meta %>% filter(hit, phenotype != "Height_ALL") %>%
		group_by(case_control) %>% 
		filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
		summarise(count = length(unique(paste(Region, phenotype))))

	# Ensure that the phenotypes going through to the plot have been meta-analysed
	meta_analysed_traits <- paste0(gsub(".*/([A-Za-z0-9]+)_.*gz", "\\1", files), "_",
		gsub(".*/[A-Za-z0-9]+_([A-Z]+).*gz", "\\1", files))

	# Damaging, unique (gene, phenotype) pairs, split by dataset and ancestry.
	dt_gene_hits_unique <- biobank_hits %>% 
		filter(paste(phenotype, sex, sep="_") %in% meta_analysed_traits) %>%
		filter(Group %in% c(
			"pLoF",
			"pLoF;damaging_missense_or_protein_altering",
			"damaging_missense_or_protein_altering")
		) %>% group_by(dataset, ancestry, case_control) %>% 
		filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

	dt_gene_hits_unique_no_height <- dt_gene_hits_unique %>% 
		filter(phenotype != "Height")%>% 
		summarise(count = length(unique(paste(Region, phenotype))))
	dt_gene_hits_unique <- dt_gene_hits_unique %>%
		summarise(count = length(unique(paste(Region, phenotype))))

	plot_unique <- rbind(
		dt_gene_hits_unique %>% mutate(dataset = unlist(renaming_plot_biobank_list[dataset])),
		meta_list_unique %>% mutate(ancestry = "Meta", dataset = "Meta")
		)
	plot_unique$dataset <- factor(plot_unique$dataset,
		levels = c(sort(setdiff(unique(plot_unique$dataset), "Meta")), "Meta"))

	plot_unique_no_height <- rbind(
		dt_gene_hits_unique_no_height %>% mutate(dataset = unlist(renaming_plot_biobank_list[dataset])),
		meta_list_unique_no_height %>% mutate(ancestry = "Meta", dataset = "Meta")
		)
	plot_unique_no_height$dataset <- factor(plot_unique_no_height$dataset,
		levels = c(sort(setdiff(unique(plot_unique_no_height$dataset), "Meta")), "Meta"))

	return(list(including_height = plot_unique, excluding_height = plot_unique_no_height))
}

tables_for_plots <- extract_hit_counts_for_figure(dt_gene_hits_all, meta_list)
tables_for_plots_cauchy <- extract_hit_counts_for_figure(dt_gene_cauchy_hits_all, meta_cauchy)

fwrite(tables_for_plots$including_height, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_data.tsv.gz", sep='\t', quote=FALSE)
fwrite(tables_for_plots$excluding_height, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_no_height_data.tsv.gz", sep='\t', quote=FALSE)

fwrite(tables_for_plots_cauchy$including_height, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_cauchy_hits_data.tsv.gz", sep='\t', quote=FALSE)
fwrite(tables_for_plots_cauchy$excluding_height, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_cauchy_hits_no_height_data.tsv.gz", sep='\t', quote=FALSE)

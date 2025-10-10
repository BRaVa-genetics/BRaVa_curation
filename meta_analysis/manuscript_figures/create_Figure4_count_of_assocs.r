#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)

source("../meta_analysis_utils.r")

# Count the number of hits, split by all of the data sets and ancestries
files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks", full.names=TRUE, recursive=TRUE)
files <- grep("/cleaned/", files, value=TRUE)
gene_files <- grep("\\.gene\\.", files, value=TRUE)

file_info <- rbindlist(lapply(gene_files, extract_file_info), fill=TRUE)
file_info <- file_info %>% mutate(filename = gene_files, dataset = gsub(".*/", "", dataset))

extract_hits <- function(filename, P_SKAT_cutoff=2.5e-7, P_SKAT_O_cutoff=2.5e-7, P_Burden_cutoff=6.7e-7) {
	dt <- fread(filename) %>% filter(
		(Pvalue < P_SKAT_O_cutoff) | 
		(Pvalue_SKAT < P_SKAT_cutoff) |
		(Pvalue_Burden < P_Burden_cutoff))
	melted_dt <- rbind(
		dt %>% select(-c("Pvalue_SKAT", "Pvalue_Burden")) %>% mutate(Test = "Pvalue") %>% filter(Pvalue < P_SKAT_O_cutoff),
		dt %>% select(-c("Pvalue", "Pvalue_Burden")) %>% rename(Pvalue = Pvalue_SKAT) %>% mutate(Test = "Pvalue_SKAT") %>% filter(Pvalue < P_SKAT_cutoff),
		dt %>% select(-c("Pvalue_SKAT", "Pvalue")) %>% rename(Pvalue = Pvalue_Burden) %>% mutate(Test = "Pvalue_Burden") %>% filter(Pvalue < P_Burden_cutoff)
		)
	return(melted_dt)
}

dt_gene_hits_all <- list()
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
	dt_gene_hits_all[[i]] <- melted_dt
}

dt_gene_hits_all_list <- dt_gene_hits_all
dt_gene_hits_all <- rbindlist(dt_gene_hits_all_list, fill=TRUE)
dt_gene_hits_all <- dt_gene_hits_all %>% filter(max_MAF %in% c("0.001", "1e-04"))
setkeyv(dt_gene_hits_all, c("phenotype", "dataset", "ancestry", "sex"))

# Code to generate inflation_summaries.tsv.gz is 'extract_genomic_inflation.r' in the folder above
# Run the code on everything using Rscript extract_genomic_inflation.r to update it.
# Here, we must remove any files that have been deemed to be inflated.
dt_inflation <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz")
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

# Remove those (phenotype, biobank, sex) files from the meta-analysis
setkeyv(dt_inflation, c("phenotype", "dataset", "ancestry", "sex"))
dt_gene_hits_all <- setdiff(dt_gene_hits_all, merge(dt_gene_hits_all, dt_inflation))

dt_gene_hits_all <- dt_gene_hits_all %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

# This is the meta-analysis results
files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100", full.names=TRUE)
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
		(Pvalue < 6.7e-7 & class == "Burden" & type == "Inverse variance weighted") |
		(Pvalue < 2.5e-7 & class %in% c("SKAT", "SKAT-O") & type == "Stouffer")))
}

meta_list <- rbindlist(meta_list) %>% mutate(case_control =
	ifelse(gsub("_.*", "", phenotype) %in% case_ctrl, TRUE,
		ifelse(gsub("_.*", "", phenotype) %in% cts, FALSE, NA)))
# Write the results
fwrite(file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_full_meta_101025.tsv.gz",
	meta_list %>% filter(hit), sep='\t', quote=FALSE)

meta_list_unique <- meta_list %>% filter(hit) %>% group_by(case_control) %>% 
	filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
	summarise(count = length(unique(paste(Region, phenotype))))

# Ensure that the phenotypes going through to the plot have been meta-analysed/
dt_gene_hits <- dt_gene_hits_all  %>% 
	filter(paste(phenotype, sex, sep="_") %in% unique(meta_list$phenotype)) %>% 
	group_by(max_MAF, Group, Test, dataset, ancestry) %>% summarise(count = n())

# Damaging, unique (gene, phenotype) pairs, split by dataset and ancestry.
dt_gene_hits_unique <- dt_gene_hits_all %>% 
	filter(paste(phenotype, sex, sep="_") %in% unique(meta_list$phenotype)) %>%
	filter(Group %in% c(
		"pLoF",
		"pLoF;damaging_missense_or_protein_altering",
		"damaging_missense_or_protein_altering")
	) %>% group_by(dataset, ancestry, case_control) %>% 
	filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>% 
	summarise(count = length(unique(paste(Region, phenotype))))

plot_unique <- rbind(
	dt_gene_hits_unique %>% mutate(dataset = unlist(renaming_plot_biobank_list[dataset])),
	meta_list_unique %>% mutate(ancestry = "Meta", dataset = "Meta")
	)
plot_unique$dataset <- factor(plot_unique$dataset,
	levels = c(sort(setdiff(unique(plot_unique$dataset), "Meta")), "Meta"))
fwrite(plot_unique, file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_data.tsv.gz", sep='\t', quote=FALSE)

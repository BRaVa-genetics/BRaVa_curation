library(data.table)
library(dplyr)

source("../meta_analysis_utils.r")

# Next, evaluate the Beta_Burden sumstats including and excluding Europe and
# compare the effect sizes

# First, let's get this working on the cluster - everything here up to reading
# the data back in is carried out on the cluster
# Read in all of the information

files <- dir(
	path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/",
	full.names=TRUE)
gene_files <-  files[!grepl("extra_cauchy", files)]

# First use the ALL or EUR file to define the set of pairs that make it through
# (this'll keep the file size down).
gene_files_all <- gene_files[grepl(".*cutoff.tsv.gz", gene_files)]

gene_files_all_list <- list()
for (file in gene_files_all)
{
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	cat(paste0(phe, "\n"))
	gene_files_all_list[[phe]] <- fread(file) %>% 
		filter(
			class == "Burden",
			type == "Inverse variance weighted", Pvalue < 6.7e-7) %>%
		mutate(phenotype = phe, ancestry = "ALL")
	setkeyv(gene_files_all_list[[phe]],
		c("Region", "Group", "max_MAF", "phenotype"))
}

gene_files <- dir(
	path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	full.names=TRUE, recursive=TRUE)
gene_files_pop <- grep("/(AFR|AMR|EAS|EUR|SAS|non_EUR)/", gene_files, value=TRUE)
gene_files_pop <- gene_files_pop[
	grep(".*cutoff.(AFR|AMR|EAS|EUR|SAS|non_EUR).tsv.gz", gene_files_pop)]

# Loop over the ancestries and loop over the phenotypes, merging at each step.
gene_files_list <- list()
for (file in gene_files_pop) {
	# Extract the phenotype name
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	ancestry <- ifelse(
		grepl(".*cutoff.([A-Za-z_]*).tsv.gz", file),
		gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file), "ALL")
	cat(paste0(phe, ": ", ancestry, "\n"))
	if (is.null(gene_files_list[[phe]])) {
		gene_files_list[[phe]] <- list()
	}
	gene_files_list[[phe]][[ancestry]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted") %>%
		mutate(phenotype = phe,
			ancestry = gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)) %>%
		rename(
			Pvalue_pop = Pvalue,
			BETA_Burden_pop = BETA_Burden,
			chisq_het_pop = chisq_het,
			Pvalue_het_pop = Pvalue_het,
			SE_Burden_pop = SE_Burden,
			BETA_meta_pop = BETA_meta,
			sum_weights_pop = sum_weights,
			ancestry_pop = ancestry,
			df_pop = df 
			) %>% select(-c("Stat", "type", "class")) 
	setkeyv(gene_files_list[[phe]][[ancestry]],
		c("Region", "Group", "max_MAF", "phenotype"))
	gene_files_list[[phe]][[ancestry]] <- merge(
		gene_files_list[[phe]][[ancestry]], gene_files_all_list[[phe]] %>% 
		select(-c("Stat", "type", "class")))
}

for (phe in names(gene_files_list)) {
	gene_files_list[[phe]] <- rbindlist(gene_files_list[[phe]])
}
gene_files_list <- rbindlist(gene_files_list)
gene_files_list <- gene_files_list %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

fwrite(gene_files_list,
	file=paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
		"Burden_comparison_for_plotting_vs_ALL.tsv.gz"), sep='\t')

# The same thing, but comparing to EUR meta

# First use the ALL or EUR file to define the set of pairs that make it through
# (this'll keep the file size down).

gene_files <- dir(
	path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	full.names=TRUE, recursive=TRUE)
gene_files_EUR <- grep("/EUR/", gene_files, value=TRUE)
gene_files_EUR <- gene_files_EUR[grep(".*cutoff.EUR.tsv.gz", gene_files_EUR)]

gene_files_EUR_list <- list()
for (file in gene_files_EUR)
{
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	cat(paste0(phe, "\n"))
	gene_files_EUR_list[[phe]] <- fread(file) %>% 
		filter(
			class == "Burden",
			type == "Inverse variance weighted",
			Pvalue < 6.7e-7) %>%
		mutate(phenotype = phe, ancestry = "EUR")
	setkeyv(gene_files_EUR_list[[phe]],
		c("Region", "Group", "max_MAF", "phenotype"))
}

gene_files_pop <- grep("/(AFR|AMR|EAS|SAS|non_EUR)/", gene_files, value=TRUE)

# Loop over the ancestries and loop over the phenotypes, merging at each step.
gene_files_list <- list()
for (file in gene_files_pop) {
	# Extract the phenotype name
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	ancestry <- gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)

	cat(paste0(phe, ": ", ancestry, "\n"))
	if (is.null(gene_files_list[[phe]])) {
		gene_files_list[[phe]] <- list()
	}

	gene_files_list[[phe]][[ancestry]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted") %>%
		mutate(phenotype = phe,
			ancestry = gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)) %>%
		rename(
			Pvalue_pop = Pvalue,
			BETA_Burden_pop = BETA_Burden,
			chisq_het_pop = chisq_het,
			Pvalue_het_pop = Pvalue_het,
			SE_Burden_pop = SE_Burden,
			BETA_meta_pop = BETA_meta,
			sum_weights_pop = sum_weights,
			ancestry_pop = ancestry,
			df_pop = df 
			) %>% select(-c("Stat", "type", "class")) 
	setkeyv(gene_files_list[[phe]][[ancestry]],
		c("Region", "Group", "max_MAF", "phenotype"))
	if (phe %in% names(gene_files_EUR_list)) {
		gene_files_list[[phe]][[ancestry]] <- merge(
			gene_files_list[[phe]][[ancestry]],
			gene_files_EUR_list[[phe]] %>% select(-c("Stat", "type", "class")))
	} else {
		gene_files_list[[phe]] <- NULL
	}
}

for (phe in names(gene_files_list)) {
	gene_files_list[[phe]] <- rbindlist(gene_files_list[[phe]])
}
gene_files_list <- rbindlist(gene_files_list)
gene_files_list <- gene_files_list %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

fwrite(gene_files_list,
	file=paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
		"Burden_comparison_for_plotting_vs_EUR.tsv.gz"),
	sep='\t')

#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)

source("../meta_analysis_utils.r")

files_list <- c("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/just_uk-biobank_and_all-of-us",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/AMR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/AFR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/EAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/EUR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/SAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/non_EUR")

meta_list <- list()
for (file in files_list) {
	files <- dir(file, full.names=TRUE)
	files <- grep(".vcf.gz$", files, value=TRUE)
	print(files)
	meta_list[[file]] <- list()

	for (f in files) {
		cat(f, "\n")
		cmd <- paste(
		"bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", f)
		dt <- fread(cmd = cmd)
		if (nrow(dt) == 0) {
			next 
		} else {
			dt <- dt %>% 
			rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7,
				`P-value`=V8) %>%
			mutate(BETA=as.numeric(BETA), SE=as.numeric(SE),
				`P-value`=-as.numeric(`P-value`))
		}
		dt <- data.table(dt) %>% filter(!is.na(`P-value`))
		setkey(dt, "P-value")
		# Filter dt remove variants on a sliding scale if the P-value is > 0.01
		dt <- dt %>% mutate(remove = ifelse(`P-value` < -2, FALSE, runif(n()) < 10^(`P-value`/2)))
		dt <- dt %>% filter(!remove) %>% select(-remove)
		meta_list[[file]][[f]] <- dt
		phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", f)
		type <- ifelse(phenotype %in% phenotype_class$continuous, "continuous", "binary")
		meta_list[[file]][[f]]$phenotype <- phenotype
		meta_list[[file]][[f]]$type <- type
	}
	meta_list[[file]] <- rbindlist(meta_list[[file]])
	meta_list[[file]]$ancestry <- ifelse(file == "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100",
		"all", gsub(".*/", "", file))
}

meta <- rbindlist(meta_list)
fwrite(meta,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/manhattan_plots_null_results_downsampled.tsv.gz",
	sep="\t", quote=FALSE)

# The file will contain *all* results for the meta-analysis
meta_files <- grep(".vcf.gz$",
	dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100",
		full.names=TRUE), value=TRUE)
dt_plot_list <- list()
for (file in meta_files)
{
	cat(file, "\n")
	cmd <- paste(
		"bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", file)
	dt <- fread(cmd = cmd) %>% 
		rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7,
			`P-value`=V8) %>%
		mutate(
			BETA=as.numeric(BETA), SE=as.numeric(SE),
			`P-value`=-as.numeric(`P-value`))
	dt <- data.table(dt) %>% filter(!is.na(`P-value`))
	setkey(dt, "P-value")
	dt_plot_list[[file]] <- dt
	phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", file)
	type <- ifelse(phenotype %in% phenotype_class$continuous, "continuous", "binary")
	dt_plot_list[[file]]$phenotype <- phenotype
	dt_plot_list[[file]]$type <- type
}

dt_manhattan_plots <- rbindlist(dt_plot_list)
fwrite(dt_manhattan_plots,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/manhattan_plots.tsv.gz",
	sep="\t", quote=FALSE)

# The code below exports all of the gene-level results in a single file. Note that this is for 
# comparing variant results against gene-level results to see when a single variant is driving
# a given signal

meta_files <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100"

files <- dir(meta_files, full.names=TRUE)
files <- grep("tsv.gz", files, value=TRUE)
print(files)
meta_list <- list()

for (f in files) {
	cat(f, "\n")
	phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", f)
	meta_list[[f]] <- fread(f) %>% filter(
		max_MAF %in% c("1e-04", "0.001"),
		Group %in% c(
			"damaging_missense_or_protein_altering",
			"pLoF",
			"pLoF;damaging_missense_or_protein_altering")) %>%
	mutate(phenotype = phenotype) %>% filter(
		(class == "Burden" & type == "Inverse variance weighted") |
		(class %in% c("SKAT", "SKAT-O") & type == "Stouffer"))
}
meta <- rbindlist(meta_list)
meta$ancestry <- "all"

fwrite(meta, sep = "\t", quote=FALSE,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/all_assocs_from_meta_051125.tsv.gz")

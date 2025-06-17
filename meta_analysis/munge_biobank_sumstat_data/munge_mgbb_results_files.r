#!/bin/Rscript
library(data.table)
library(dplyr)

biobank <- "mgbb"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- FALSE

cloud_gene_folder <- "20240515"
cloud_variant_folder <- "20240515"

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
		cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	)
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
		cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
	)
}

source("../meta_analysis_utils.r")

rename_files <- function(results_type="gene", analysis_name="pilot")
{
	# Change the filenames with gsub
	from <- list.files(paste0(data_dir, "/", results_type), pattern=analysis_name, full.names=TRUE)
	start <- gsub(paste0("(.*", analysis_name, ".*)\\.(ALL|F|M).*"), "\\1", from)
	end <- gsub(paste0("(.*", analysis_name, ".*)(\\.(ALL|F|M).*$)"), "\\2", from)
	for (i in 1:length(from)) {
		if (!(extract_file_info(from[i])$sex %in% c("ALL", "M", "F"))) {
			file.rename(from[i], paste0(start[i], ".1", end[i]))
		}
	}
}

rename_files()
rename_files(results_type="variant")

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/variant",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)

# Continous phenotypes provided in September

cloud_gene_folder <- "20240913"
cloud_variant_folder <- "20240913"

system(paste0("mkdir -p ", data_dir, "/gene/", cloud_gene_folder))
system(paste0("mkdir -p ", data_dir, "/variant/", cloud_variant_folder))

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
		cloud_gene_folder, "/*_Median.*.gene.*gz ", data_dir, "/gene/", cloud_gene_folder)
	)
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
		cloud_variant_folder, "/*_Median.*.variant.*gz ", data_dir, "/variant/", cloud_variant_folder)
	)
}

# Note that a new function is required for the renaming here
rename_files <- function(results_type="gene", analysis_name="pilot", subfolder="20240913")
{
	# Change the filenames with gsub
	from <- list.files(paste0(data_dir, "/", results_type, "/", subfolder),
		pattern=analysis_name, full.names=TRUE)
	start <- gsub("_Median", "", gsub(paste0("(.*", analysis_name, ".*)\\.(ALL|F|M).*"), "\\1", from))
	end <- gsub(paste0("(.*", analysis_name, ".*)(\\.(ALL|F|M).*$)"), "\\2", from)
	# These are the continuous files - so there should only be a single N
	end <- gsub("(.*)\\.([0-9]+)\\.([0-9]+)\\.(.*txt.gz)", "\\1.\\3.\\4", end) 
	for (i in 1:length(from)) {
		if (!(extract_file_info(from[i])$sex %in% c("ALL", "M", "F"))) {
			file.rename(from[i], paste0(start[i], ".1", end[i]))
		}
	}
}

rename_files()
rename_files(results_type="variant")

for (file in dir(paste0(data_dir, "/gene/", cloud_gene_folder), full.names=TRUE)) {
	dt <- fread(file, nrows=0)
	if (all(names(dt) == paste0("V", seq(1, 11)))) {
		dt <- fread(file)
		names(dt) <- c("Region", "Group", "max_MAF",
			"Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
			"BETA_Burden", "SE_Burden", "MAC", "Number_rare",
			"Number_ultra_rare")
		cat("Writing file, pasting in column names\n")
		fwrite(dt, sep='\t', quote=FALSE, file=file)
	}
}

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/gene/20240913",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/variant/20240913",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)

#!/bin/Rscript
library(data.table)
library(dplyr)

biobank <- "biome"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

cloud_gene_folder <- "gene_run1"
cloud_variant_folder <- "variant"

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
		file_info <- as.list(strsplit(from[i], split="\\.")[[1]])
		if (!(file_info[[6]] %in% c("ALL", "M", "F"))) {
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

#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "all-of-us"

data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")
system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- FALSE

source("meta_analysis_utils.r")

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/v8/*.gene.* ",
		data_dir, "/gene/")
	)
}

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

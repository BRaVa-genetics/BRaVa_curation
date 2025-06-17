#!/bin/Rscript
library(data.table)
library(dplyr)

biobank <- "pmbb"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
		"/RUN2/*.gene.* ", data_dir, "/gene/")
	)
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
		"/RUN2/*.variant.* ", data_dir, "/variant/")
	)
}

# Remove the MatHem file that displays inflation
system(paste0("rm ", data_dir, "/gene/pmbb.rodriguez.pilot.MatHem.01.FEMALE.EUR.589.12845.SAIGE.gene.20240515.txt.gz"))
system(paste0("rm ", data_dir, "/variant/pmbb.rodriguez.pilot.MatHem.01.FEMALE.EUR.589.12845.SAIGE.variant.20240515.txt.gz"))
source("../meta_analysis_utils.r")

rename_files <- function(results_type="gene", analysis_name="pilot")
{
	# Change the filenames with gsub
	from <- list.files(paste0(data_dir, "/", results_type), pattern=analysis_name, full.names=TRUE)
	start <- gsub(paste0("(.*", analysis_name, ".*)\\.(ALL|F|M).*"), "\\1", from)
	end <- gsub(paste0("(.*", analysis_name, ".*)(\\.(ALL|F|M).*$)"), "\\2", from)
	for (i in 1:length(from)) {
		if (!(extract_file_info(from[i])$sex %in% c("ALL", "M", "F", "all", "male", "female", "FEMALE", "MALE"))) {
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

#!/bin/Rscript

biobank <- "genes-and-health"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

freeze_to_match <- "JULY23Freeze"
round_to_match <- "Round1"

rename_files <- function(results_type="gene")
{
	# Change the filenames with gsub
	from <- list.files(paste0(data_dir, "/", results_type), pattern=freeze_to_match, full.names=TRUE)
	from <- from[which(grepl(paste0(round_to_match, ".*", freeze_to_match), from))]
	phenotype <- gsub(paste0(".*", round_to_match, "\\.(.*)\\.", freeze_to_match, ".*"), "\\1", from)
	new_phenotype <- gsub("\\.", "_", phenotype)
	new_phenotype <- gsub("[_]+", "_", new_phenotype)
	new_phenotype <- gsub("_$", "", new_phenotype)
	start <- gsub(paste0("(.*", round_to_match, "\\.).*"), "\\1", from)
	end <- gsub(paste0(".*(\\.", freeze_to_match, ".*)"), "\\1", from)
	cat(paste0("from: ", from, "\nto: ", paste0(start, new_phenotype, end)))
	file.rename(from, paste0(start, new_phenotype, end))
}

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/*.gene.* ",
		data_dir, "/gene/")
	)
	system(paste0("rm ", data_dir, "/gene/PRELIMINARY_DATA_NOT_FINAL_GNH_SUBMISSION_*"))
	rename_files()
	
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/*.variant.* ",
		data_dir, "/variant/")
	)
	system(paste0("rm ", data_dir, "/variant/PRELIMINARY_DATA_NOT_FINAL_GNH_SUBMISSION_*"))
	rename_files(results_type="variant")
}

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

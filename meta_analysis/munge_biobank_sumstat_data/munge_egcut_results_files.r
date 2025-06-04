#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "egcut"
last_name <- "Joeloo"
analysis_name <- "pilot"
freeze_number <- "2024v03"

data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

cloud_gene_folder <- NULL
cloud_variant_folder <- NULL

if (download) {
	gene_cmd_string <- ifelse(
		is.null(cloud_gene_folder),
		paste0("gsutil -m cp gs://brava-meta-upload-",
			biobank, "/*.gene.* ", data_dir, "/gene/"),
		paste0("gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")

		)
	system(gene_cmd_string)
	
	# variant_cmd_string <- ifelse(
	# 	is.null(cloud_variant_folder),
	# 	paste0("gsutil -m cp gs://brava-meta-upload-",
	# 		biobank, "/*.variant.* ", data_dir, "/variant/"),
	# 	paste0("gsutil -m cp gs://brava-meta-upload-", biobank,
	# 		"/", cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")

	# 	)
	# system(variant_cmd_string)
}

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

# system(paste0("Rscript munge_results_files_Group_names.r",
# 	" --folder ", data_dir, "/variant",
# 	" --type ", "variant",
# 	" --write",
# 	" --out_folder ", out_data_dir, "/variant")
# )

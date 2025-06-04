#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "gel"
last_name <- "kousathanas"
analysis_name <- "pilot"
freeze_number <- "1"

data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

# cloud_gene_folder <- "pilot_32_phenotypes_20240508"
# cloud_variant_folder <- "pilot_32_phenotypes_20240508"

# if (download) {
# 	system(paste0(
# 		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
# 		cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
# 	)
# 	system(paste0(
# 		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
# 		cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
# 	)
# }

# system(paste0("Rscript munge_results_files_Group_names.r",
# 	" --folder ", data_dir, "/gene",
# 	" --type ", "gene",
# 	" --write",
# 	" --out_folder ", out_data_dir, "/gene")
# )

# system(paste0("Rscript munge_results_files_Group_names.r",
# 	" --folder ", data_dir, "/variant",
# 	" --type ", "variant",
# 	" --write",
# 	" --out_folder ", out_data_dir, "/variant")
# )

cloud_gene_folder <- "pilot_32_phenotypes_allancestries_20250217/RVAT_summaries/"
cloud_gene_folders <- paste0(cloud_gene_folder, c("AFR", "EAS", "EUR", "SAS"))
cloud_variant_folders <- cloud_gene_folders

for (cloud_gene_folder in cloud_gene_folders) {
	if (download) {
		system(paste0(
			"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
			cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
		)
	}
}

# for (cloud_variant_folder in cloud_variant_folders) {
# 	system(paste0(
# 			"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
# 			cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
# 		)
# }

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

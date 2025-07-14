#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "gel"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

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

for (cloud_gene_folder in cloud_gene_folders) {
	if (download) {
		system(paste0(
			"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
			cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
		)
	}
}

cloud_variant_folder <- "pilot_32_phenotypes_allancestries_20250217/RVAT_single_variant_summaries/"
cloud_variant_folders <- paste0(cloud_variant_folder, c("AFR", "EAS", "EUR", "SAS"))

for (cloud_variant_folder in cloud_variant_folders) {
	if (download) {
		system(paste0(
				"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
				cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
			)
	}
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

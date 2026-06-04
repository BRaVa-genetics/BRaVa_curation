#!/bin/Rscript
library(data.table)
library(dplyr)

biobank <- "bbj"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

cloud_gene_folder <- "pilot_44phenotypes/binary_traits"
cloud_variant_folder <- "pilot_44phenotypes/binary_traits"

if (download) {
	# system(paste0(
	# 	"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
	# 	cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	# )
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/",
		cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
	)
}

# system(paste0("Rscript munge_results_files_Group_names.r",
# 	" --folder ", data_dir, "/gene",
# 	" --type ", "gene",
# 	" --write",
# 	" --out_folder ", out_data_dir, "/gene")
# )

system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", data_dir, "/variant",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)

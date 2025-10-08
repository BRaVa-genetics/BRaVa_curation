#!/bin/Rscript
library(data.table)
library(dplyr)

biobank <- "ccpm"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

cloud_gene_folder <- "brava-results/august_2025_CCPM_results/binary_gene_based"
cloud_variant_folder <- "brava-results/august_2025_CCPM_results/binary_variant_based"

# As a result of the above, the remainder here is hard-coded
if (download) {
	gene_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	system(gene_cmd_string)
	
	variant_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
	system(variant_cmd_string)
}

# Remove the MatHem file that displays inflation
system(paste0("rm ", data_dir, "/gene/CCPM.White.pilot.MatHem.exf3.FEMALE.EUR.*.txt.gz"))
system(paste0("rm ", data_dir, "/variant/CCPM.White.pilot.MatHem.exf3.FEMALE.EUR.*.txt.gz"))

# Subsequently, the continuous traits were added
cloud_gene_folder <- "brava-results/august_2025_CCPM_results/continuous_gene_based"
cloud_variant_folder <- "brava-results/august_2025_CCPM_results/continuous_variant_based"

if (download) {
	gene_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	system(gene_cmd_string)

	variant_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_variant_folder, "/*.variant.* ", data_dir, "/variant/")
	system(variant_cmd_string)
}

# Remove the 'MIXED' analyses
system(paste0("rm ", data_dir, "/gene/CCPM.*.MIXED.*"))
system(paste0("rm ", data_dir, "/variant/CCPM.*.MIXED.*"))

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

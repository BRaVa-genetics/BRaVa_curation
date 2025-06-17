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

# DEV: need to fix this - their gene level results are variant level results
cloud_gene_folder <- "brava-results/gene_based_feb_2025_16"
# The following is odd because these files are the variant outputs from a SAIGE-gene call
cloud_variant_folder <- "brava-results/gene_based"

# As a result of the above, the remainder here is hard-coded
if (download) {
	gene_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	system(gene_cmd_string)
	
	variant_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_variant_folder, "/*.gene.* ", data_dir, "/variant/")
	system(variant_cmd_string)
}

# Remove the MatHem file that displays inflation
system(paste0("rm ", data_dir, "/gene/CCPM.White.pilot.MatHem.f4.FEMALE.EUR.1925.29807.SAIGE.gene.20250216.txt.gz"))
system(paste0("rm ", data_dir, "/variant/CCPM.White.pilot.MatHem.f4.FEMALE.EUR.1925.29807.SAIGE.gene.20250210.txt.gz"))

# Rename the variant files.
for (file in dir(paste0(data_dir, "/variant"), full.names=TRUE)) {
	system(paste("mv", file, gsub("gene.20250210.txt.gz", "variant.20250210.txt.gz", file)))
}

# Subsequently, the continuous traits were added
cloud_gene_folder <- "brava-results/gene_based_continuous"

if (download) {
	gene_cmd_string <- paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank,
			"/", cloud_gene_folder, "/*.gene.* ", data_dir, "/gene/")
	system(gene_cmd_string)
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



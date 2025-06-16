#!/bin/Rscript
library(data.table)
library(dplyr)

# For naming files
biobank <- "all-of-us"
last_name <- "Lu"
analysis_name <- "pilot"
freeze_number <- "JULY23Freeze"

data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")
system(paste0("mkdir -p ", data_dir, "/gene"))
system(paste0("mkdir -p ", data_dir, "/variant"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- FALSE

# DEV: fix this so that it is generated directly from the files uploaded rather than relying on a hard-coded file from someone else.
data_dictionary <- fread("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/brava_44_pilot_phenotype_filename.csv")
# Add the AST files that were missed
# Note that CRP is a binary trait in all of us, so is not included
data_dictionary <- rbind(data_dictionary,
	data.table(
		"aou_code" = rep("3013721", 4),
		"description" = rep("Aspartate aminotransferase", 4),
		"brava_code" = rep("AST", 4),
		"label" = rep("AST_3013721", 4),
		"pop" = c("afr", "amr", "eur", "sas"),
		"pheno_sex" = rep("ALL", 4),
		"n_cases" = c(22276, 20472, 67636, 1299),
		"n_controls" = rep(NA, 4),
		"trait_type" = rep("continuous", 4),
		"gene_file_name" = c("all-of-us.Lu.pilot.AST_3013721.ALL.AFR.22276.SAIGE.gene.20240424.txt.gz",
							 "all-of-us.Lu.pilot.AST_3013721.ALL.AMR.20472.SAIGE.gene.20240424.txt.gz",
							 "all-of-us.Lu.pilot.AST_3013721.ALL.EUR.67636.SAIGE.gene.20240424.txt.gz",
							 "all-of-us.Lu.pilot.AST_3013721.ALL.SAS.1299.SAIGE.gene.20240424.txt.gz"),
		"var_file_name" = c("all-of-us.Lu.pilot.AST_3013721.ALL.AFR.22276.SAIGE.variant.20240424.txt.gz",
							"all-of-us.Lu.pilot.AST_3013721.ALL.AMR.20472.SAIGE.variant.20240424.txt.gz",
							"all-of-us.Lu.pilot.AST_3013721.ALL.EUR.67636.SAIGE.variant.20240424.txt.gz",
							"all-of-us.Lu.pilot.AST_3013721.ALL.SAS.1299.SAIGE.variant.20240424.txt.gz")
		), fill=TRUE
	)

source("meta_analysis_utils.r")

# # These are the files that would have been retained if we just picked the trait with the largest sample size
# files_to_retain_old <- data.table(data_dictionary %>% 
# 	group_by(brava_code, pheno_sex, pop) %>% 
# 	filter(n_sample_label == max(n_sample_label)))

# These are the files that we retain using phecode x where available, otherwise phecode.
# Combine with the continuous traits
files_to_retain <- rbind(
	data.table(data_dictionary %>% filter(aou_code %in% unlist(aou_curated))),
	data.table(data_dictionary %>% filter(trait_type == "continuous"))
	)

# Note that the following files are the same (as ascertained by md5sum):
# The reason for this is likely that Wenhan uploaded some files and then after the
# fact, noticed that our requested phenotype labels didn't quite match.

# Total cholesterol
# all-of-us.Lu.pilot.cholesterol.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.TChol_3027114.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# HDL cholesterol
# all-of-us.Lu.pilot.HDL.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.HDLC_3007070.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# LDL cholesterol
# all-of-us.Lu.pilot.LDL.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz
# all-of-us.Lu.pilot.LDLC_3028288.ALL.{AFR,AMR,EAS,EUR,SAS}.{count}.SAIGE.gene.20240401.txt.gz

# T2D phenotypes
# all-of-us.Lu.pilot.T2Diab_EM_202.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2D_EM_202.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2Diab_250.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz
# all-of-us.Lu.pilot.T2D_250.2.ALL.{AFR,AMR,EAS,EUR,SAS}.{case_count}.{control_count}.SAIGE.gene.20240424.txt.gz

if (download) {
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/v7/*.gene.* ",
		data_dir, "/gene/")
	)
	system(paste0(
		"gsutil -m cp gs://brava-meta-upload-", biobank, "/v7/*.variant.* ",
		data_dir, "/variant/")
	)
}

files_to_retain <- files_to_retain %>% 
	mutate(new_gene_file_name = paste(
		biobank,
		last_name,
		analysis_name,
		brava_code,
		freeze_number,
		pheno_sex,
		toupper(pop),
		ifelse(trait_type=="continuous", n_cases, paste(n_cases, n_controls, sep=".")),
		"SAIGE",
		"gene",
		gsub("^.*gene.([0-9]+.*$)", "\\1", gene_file_name),
		sep=".")) %>% 
	mutate(new_variant_file_name = gsub(
		"\\.gene\\.", "\\.variant\\.", new_gene_file_name))

from <- paste0(data_dir, "/gene/", files_to_retain$gene_file_name)
to <- paste0(out_data_dir, "/gene/", files_to_retain$new_gene_file_name)

file.rename(from, to)
system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", out_data_dir, "/gene",
	" --type ", "gene",
	" --write",
	" --out_folder ", out_data_dir, "/gene")
)

from <- paste0(data_dir, "/variant/", files_to_retain$var_file_name)
to <- paste0(out_data_dir, "/variant/", files_to_retain$new_variant_file_name)

file.rename(from, to)
system(paste0("Rscript munge_results_files_Group_names.r",
	" --folder ", out_data_dir, "/variant",
	" --type ", "variant",
	" --write",
	" --out_folder ", out_data_dir, "/variant")
)


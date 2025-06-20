#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("meta_analysis_utils.r")
source("../phenotypes/BRaVa_phenotypes_utils.r")

main <- function(args)
{
	data_dir <- args$vcf_data_dir
	n_cases <- args$n_cases
	out_dir <- args$out_dir
	phe <- args$phenotypeID

	out_meta_results_dir <- paste0(out_dir, "/variant/n_cases_", n_cases)
	
	# Ensure that the folder is already present
	system(paste("mkdir -p", out_meta_results_dir))

	# Assumes that we have the files locally within a file structure as defined
	# in the munging scripts.
	# First, let's determine the collection of phenotypes that we are testing
	# and the collection of (population, phenotype) pairs for that biobank to 
	# include in the meta-analysis

	biobanks <- dir(data_dir)[file.info(dir(data_dir, full.names=TRUE))$isdir]
	biobanks <- intersect(biobanks, names(file_check_information$dataset))
	results_dt_list <- list()

	vcf_variant_files <- dir(data_dir)
	vcf_variant_files <- vcf_variant_files[-grep(".gz.csi$",
		vcf_variant_files)]
	results <- lapply(vcf_variant_files, extract_file_info)
	results_dt <- data.table(
		filename = paste0(data_dir, "/", vcf_variant_files),
		phenotypeID = sapply(results, `[[`, "phenotype"),
		pop = sapply(results, `[[`, "ancestry"),
		biobank = sapply(results, `[[`, "dataset"),
		sex = sapply(results, `[[`, "sex")
	)

	if (is.null(phe)) {
		phes <- BRaVa_pilot_phenotypes
	} else {
		phes <- phe
	}

	# Here, we must remove any files that have been deemed to be inflated.
	dt_inflation <- fread(args$inflation_file)
	dt_inflation <- unique(dt_inflation %>% filter(Group == "synonymous") %>% 
		filter(max_MAF != 0.01,
			lambda_value > 1.3,
			!(lambda_type %in% c(
				"lambda_50_Burden", "lambda_50_SKAT", "lambda_50"))) %>% 
		select(phenotype, dataset, ancestry, sex))
	# Manual curation, adding the following (biobank, trait) tuples containing
	# spurious associations
	dt_inflation <- rbind(dt_inflation, data.table(
		phenotype = c("ColonRectCanc", "Height"),
		dataset = c("egcut", "mgbb"),
		ancestry = c("EUR", "AMR"),
		sex = c("ALL", "ALL"))
	) %>% rename(phenotypeID = phenotype, biobank = dataset, pop = ancestry)

	# Remove those (phenotype, biobank, sex) files from the meta-analysis
	setkeyv(dt_inflation, c("phenotypeID", "biobank", "pop", "sex"))
	setkeyv(results_dt,  c("phenotypeID", "biobank", "pop", "sex"))
	results_dt <- setdiff(results_dt, merge(dt_inflation, results_dt))

	# Everything
	for (phe in phes) {
		for (s in c("ALL", "M", "F")) {
			files_vcf <- (results_dt %>% 
				filter(phenotypeID == phe, sex == s))$filename
			# For each vcf file - ensure that it is not empty.
			empty_vcfs <- which(
				sapply(files_vcf, function(f) file.info(f)$size == 0))
			if (length(empty_vcfs) > 0) {
				files_vcf <- files_vcf[-empty_vcfs]
			}
			if (length(files_vcf) <= 1) { 
				cat(paste("Either the phenotype is not present, or there is",
					"only a single file for:\n"))
				cat(phe, s, "\n")
			} else {
				files_vcf <- paste(files_vcf, collapse=",")
				out <- paste0(out_meta_results_dir, "/", phe, "_", s,
					"_variant_meta_analysis_", n_cases, "_cutoff.vcf.gz")
				cat(paste0("carrying out meta-analysis of ",
					phe, " in ", s, "\n"))
				cat(paste0("\nFiles in the analysis: ",
					paste0(strsplit(files_vcf, split=",")[[1]], collapse='\n'),
					"\n"))
				system(paste(
					"sbatch",
					"run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh",
					files_vcf, out))
				cat(paste0("submitted meta-analysis of ", phe,
					" completed\n\n"))
			}
		}
	}

	# We also want to run superpopulation specific meta-analysis, non-EUR,
	# and leave-one-biobank-out meta-analysis as well

	# For each superpopulation
	for (phe in phes) {
		for (s in c("ALL", "M", "F")) {
			
			files_vcf <- (results_dt %>%
				filter(phenotypeID == phe, sex == s))$filename
			files_info <- lapply(files_vcf, extract_file_info)
			to_subset <- data.table(
				filename = files_vcf,
				phenotypeID = sapply(files_info, `[[`, "phenotype"),
				pop = sapply(files_info, `[[`, "ancestry"),
				sex = sapply(files_info, `[[`, "sex")
			)

			for (p in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
				files_vcf_tmp <- (to_subset %>% filter(pop == !!p))$filename
				# For each vcf file - ensure that it is not empty.
				empty_vcfs <- which(
					sapply(files_vcf_tmp, function(f) file.info(f)$size == 0))
				if (length(empty_vcfs) > 0) {
					files_vcf_tmp <- files_vcf_tmp[-empty_vcfs]
				}
				if (length(files_vcf_tmp) <= 1) { 
					cat("Either the phenotype is not present, or there is",
						"only a single file for:\n")
					cat(phe, s, p, "\n")
				} else {
					files_vcf_tmp <- paste(files_vcf_tmp, collapse=",")
					# Ensure the folder is present
					system(paste0("mkdir -p ", out_meta_results_dir, "/", p))
					out <- paste0(out_meta_results_dir, "/", p, "/",
						phe, "_", s, "_variant_meta_analysis_",
						n_cases, "_cutoff.", p, ".vcf.gz")
					cat(paste0("carrying out meta-analysis of ", phe,
						" in ", s, " for ", p, "\n"))
					cat(paste0("\nFiles in the analysis: ",
						paste0(strsplit(files_vcf_tmp, split=",")[[1]],
							collapse='\n'), "\n"))
					system(paste(
						"sbatch",
						"run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh",
						files_vcf_tmp, out))
					cat(paste0("submitted meta-analysis of ",
						phe, ":", p, " completed\n\n"))
				}
			}
		}
	}

	# for non-EUR
	for (phe in phes) {
		for (s in c("ALL", "M", "F")) {

			files_vcf <- (results_dt %>%
				filter(phenotypeID == phe, sex == s))$filename
			files_info <- lapply(files_vcf, extract_file_info)
			to_subset <- data.table(
				filename = files_vcf,
				phenotypeID = sapply(files_info, `[[`, "phenotype"),
				pop = sapply(files_info, `[[`, "ancestry"),
				sex = sapply(files_info, `[[`, "sex")
			)
			
			files_vcf_tmp <- (to_subset %>% filter(pop != "EUR"))$filename
			# For each vcf file - ensure that it is not empty.
			empty_vcfs <- which(
				sapply(files_vcf_tmp, function(f) file.info(f)$size == 0))
			if (length(empty_vcfs) > 0) {
				files_vcf_tmp <- files_vcf_tmp[-empty_vcfs]
			}
				
			if (length(files_vcf_tmp) <= 1) { 
				cat("Either the phenotype is not present, or there is",
				"only a single file for:\n")
				cat(phe,  s, "non-EUR\n")
			} else {
				files_vcf_tmp <- paste(files_vcf_tmp, collapse=",")
				# Ensure the folder is present
				system(paste0("mkdir -p ", out_meta_results_dir, "/non_EUR"))
				out <- paste0(out_meta_results_dir, "/non_EUR/",
					phe, "_", s, "_variant_meta_analysis_",
					n_cases, "_cutoff.non_EUR.vcf.gz")
				cat(paste0("carrying out meta-analysis of ", phe,
					" in ", s, " for non-EUR\n"))
				cat(paste0("\nFiles in the analysis: ",
					paste0(strsplit(files_vcf_tmp, split=",")[[1]],
						collapse='\n'), "\n"))
				system(paste(
					"sbatch",
					"run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh",
					files_vcf_tmp, out))
				cat(paste0("submitted meta-analysis of ",
					phe, ":non_EUR completed\n\n"))
			}
		}
	}

	# And now, do a leave one biobank out meta-analysis.
	# Here, we restrict 'leave one out' biobanks to those containing the
	# biobank in the original analysis
	for (b in unique(results_dt$biobank)) {
		results_dt_tmp <- results_dt %>% filter(biobank != b)
		for (phe in phes) {
			for (s in c("ALL", "M", "F")) {

				files_vcf <- (results_dt_tmp %>% 
					filter(phenotypeID == phe, sex == s))$filename
				# For each vcf file - ensure that it is not empty.
				empty_vcfs <- which(
					sapply(files_vcf, function(f) file.info(f)$size == 0))
				if (length(empty_vcfs) > 0) {
					files_vcf <- files_vcf[-empty_vcfs]
				}
				
				if (length(files_vcf) <= 1) { 
					cat("Either the phenotype is not present, or there is",
					"only a single file for:\n")
					cat(phe, s, "for all biobanks except", b, "\n")
				} else {
					files_vcf <- paste(files_vcf, collapse=",")
					# Ensure the folder is present
					system(paste0("mkdir -p ", out_meta_results_dir,
						"/minus_", b))
					out <- paste0(out_meta_results_dir, "/minus_", b, "/",
						phe, "_", s, "_variant_meta_analysis_",
						n_cases, "_cutoff.minus_", b, ".vcf.gz")
					cat(paste0("carrying out meta-analysis of ", phe,
						" in ", s, " for all biobanks except ", b, "\n"))
					cat(paste0("\nFiles in the analysis: ",
						paste0(strsplit(files_vcf, split=",")[[1]],
							collapse='\n'), "\n"))
					system(paste(
						"sbatch",
						"run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh",
						files_vcf, out))
					cat(paste0("submitted meta-analysis of ",
						phe, ":minus ", b, " completed\n\n"))
				}
			}
		}
	}

	# Finally, run a meta-analysis using just UK-biobank and All of Us
	results_dt <- results_dt %>% 
		filter(biobank %in% c("uk-biobank", "all-of-us"))
	# for non-EUR
	for (phe in phes) {
		for (s in c("ALL", "M", "F")) {

			files_vcf <- (results_dt %>% 
				filter(phenotypeID == phe, sex == s))$filename
			# For each vcf file - ensure that it is not empty.
			empty_vcfs <- which(
				sapply(files_vcf, function(f) file.info(f)$size == 0))
			if (length(empty_vcfs) > 0) {
					files_vcf <- files_vcf[-empty_vcfs]
			}
			
			if (length(files_vcf) <= 1) { 
				cat("Either the phenotype is not present, or there is",
				"only a single file for:\n")
				cat(phe,  s, "just uk-biobank and all-of-us\n")
			} else {
				files_vcf <- paste(files_vcf, collapse=",")
				# Ensure the folder is present
				system(paste0("mkdir -p ", out_meta_results_dir,
					"/just_uk-biobank_and_all-of-us"))
				out <- paste0(out_meta_results_dir,
					"/just_uk-biobank_and_all-of-us/", phe, "_", s,
					"_gene_meta_analysis_", n_cases,
					"_cutoff.just_uk-biobank_and_all-of-us.tsv.gz")
				cat(paste0("carrying out meta-analysis of ",
					phe, " in ", s, " for non-EUR\n"))
				cat(paste0("\nFiles in the analysis: ",
					paste0(strsplit(files_vcf, split=",")[[1]],
						collapse='\n'), "\n"))
				system(paste(
					"sbatch",
					"run_meta_analysis_bcftools_metal_variant_gcloud_bmrc.sh",
					files_vcf, out))
				cat(paste0("submitted meta-analysis of ",
					phe, ": just uk-biobank and all-of-us completed\n\n"))
			}
		}
	}

}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--vcf_data_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/vcf/variant",
	required=FALSE, help="Location of the folder containing the input into the meta-analysis")
parser$add_argument("--n_cases", default=100, required=FALSE,
	help="Minimum number of cases")
parser$add_argument("--out_dir", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs",
	required=FALSE, help="Output folder path")
parser$add_argument("--inflation_file", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz",
	required=FALSE, help="Inflation file")
parser$add_argument("--phenotypeID", required=FALSE, default=NULL,
	help="The phenotype ID to run meta-analysis on. Note: thus exactly must match the naming in input folder.")
args <- parser$parse_args()

main(args)

#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

source("meta_analysis_utils.r")
source("determine_AF.r")

main <- function(args)
{	
	file <- args$file_path
	out_data_dir <- args$out_data_dir
	dt <- fread(file)

	# Remove all variants for which the p.value.NA < 0.01 but the SPA has not 
	# converged
	if ("Is.SPA" %in% colnames(dt)) {
		dt <- dt %>% filter(!(SE < 1e-15 & !Is.SPA))
		dt <- data.table(dt)
	}

	if (!(any(c("N_case", "N_ctrl") %in% names(dt)))) {
		info <- extract_file_info(file)
		if ("n_cases" %in% names(info)) {
			dt$N_case <- as.integer(info$n_cases)
			dt$N_ctrl <- as.integer(info$n_controls)
		}
	}

	if (!("AF_Allele2" %in% names(dt))) {
		dt$AF_Allele2 <- NA
	} 

	if (any(is.na(dt$AF_Allele2)) & ("N_case" %in% names(dt))) {
		setDT(dt)  # dt has BETA, SE, N_case, N_ctrl
		dt <- reconstruct_df(dt)
		dt <- dt %>% mutate(
			AF_Allele2 = ifelse(!is.na(AF_Allele2), AF_Allele2, AF),
			AC_Allele2 = ifelse(!is.na(AC_Allele2), AC_Allele2, ceiling(AC))
			) %>% select(-c("AC_case", "AC_ctrl", "AC", "AF"))
	}

	setkeyv(dt, c("CHR", "POS", "Allele1", "Allele2"))
	file_tmp <- paste0("/tmp/", gsub(".txt.gz", ".tmp.txt.gz", basename(file)))
	file_vcf_tmp <- paste0("/tmp/", gsub(".txt.gz", ".vcf.gz", basename(file)))
	fwrite(dt, sep='\t', quote=FALSE, file=file_tmp)
	cmd <- paste0(
		"bcftools +munge --no-version -f /well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ",
		"-c SAIGE ", file_tmp, " -o ", file_vcf_tmp, " -Oz -W")
	cat(paste(cmd, "\n"))
	system(cmd)
	cmd <- paste0(
		"bcftools norm -f ",
		"/well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ",
		"--rm-dup both ", file_vcf_tmp, " -o ",
		gsub(".txt.gz$", ".vcf.gz", paste0(out_data_dir, "/", basename(file))), " -Oz -W"
		)
	cat(paste(cmd, "\n"))
	system(cmd)
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--file_path", default=NULL, required=TRUE,
    help="file path of file to be munged")
parser$add_argument("--out_data_dir",
	default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/vcf/variant",
	required=FALSE, help="Output file directory")
args <- parser$parse_args()

main(args)
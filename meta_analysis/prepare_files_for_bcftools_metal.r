#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)

main <- function(args)
{	
	file <- args$file_path
	out_data_dir <- args$out_data_dir
	dt <- fread(file)
	setkeyv(dt, c("CHR", "POS", "Allele1", "Allele2"))
	file_tmp <- paste0("/tmp/", gsub(".txt.gz", ".tmp.txt.gz", basename(file)))
	fwrite(dt, sep='\t', quote=FALSE, file=file_tmp)
	cmd <- paste0(
		"bcftools +munge --no-version -f /well/lindgren/dpalmer/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ",
		"-c SAIGE ", file_tmp, " -o ",
		gsub(".txt.gz$", ".vcf.gz", paste0(out_data_dir, "/", basename(file))), " -Oz -W")
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
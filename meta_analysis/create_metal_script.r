library(argparse)

main <- function(args)
{
	cat(paste0("current phenotype: ", args$phenotypeID, "\n"))
	filename <- paste0(args$out_folder, "/METAL_", args$phenotypeID, ".txt")
	writeLines(c(
		"MARKERLABEL MarkerID",
		"ALLELELABELS Allele1 Allele2",
		"EFFECTLABEL BETA",
		"STDERRLABEL SE",
		"PVALUELABEL p.value",
		"SCHEME STDERR"), filename
	)

	fileConn <- file(filename, 'a')
	# String separate the files
	for (file in strsplit(args$files, split=",")[[1]]) {
		writeLines(paste("PROCESS", file), fileConn)
	}

	writeLines(c(
		"",
		paste0("OUTFILE ", args$out_meta_analysis_folder, "/", args$phenotypeID, "_variant_meta_analysis .tbl"),
		"ANALYZE HETEROGENEITY",
		"",
		"QUIT"), fileConn
	)
	close(fileConn)
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--phenotypeID", required=TRUE,
    help="phenotypeID for the output filename")
parser$add_argument("--files", default=NULL, required=TRUE,
    help="comma separated string of files to use in the meta-analysis")
parser$add_argument("--out_folder", default=".", required=FALSE)
parser$add_argument("--out_meta_analysis_folder", default=".", required=FALSE)
args <- parser$parse_args()

main(args)

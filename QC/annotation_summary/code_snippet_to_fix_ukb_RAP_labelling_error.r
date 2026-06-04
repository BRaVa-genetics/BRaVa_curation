# Small code snippet to fix the labelling in the variant frequency annotations

library(data.table)
library(dplyr)

biobank <- "ukb_wes_450k"
count_directory <- "/Users/duncan/Repositories/BRaVa_curation/counts/"
file_grep <- paste0(biobank, ".[A-Z]{3}.chr@.BRaVa_annotations_.*_summary.tsv.gz")
files_grep <- gsub("@", "[0-9]*X?", file_grep)
files <- grep(files_grep, dir(count_directory, full.names=TRUE), value=TRUE)

for (f in files) {
	dt <- fread(f)
	dt$bin[which(dt$bin == "<0.01%")] <- "<0.1%"
	dt$bin[which(dt$bin == "0.01-0.1%")] <- "0.1-1%"
	fwrite(dt, file=f, sep="\t", quote=FALSE)
}


biobank <- "ukb_wes_450k"
count_directory <- "/Users/duncan/Repositories/BRaVa_curation/counts/"
file_grep <- paste0(biobank, ".[A-Z]{3}.chr@.vep_annotations_.*_summary.tsv.gz")
files_grep <- gsub("@", "[0-9]*X?", file_grep)
files <- grep(files_grep, dir(count_directory, full.names=TRUE), value=TRUE)

for (f in files) {
	dt <- fread(f)
	dt$bin[which(dt$bin == "<0.01%")] <- "<0.1%"
	dt$bin[which(dt$bin == "0.01-0.1%")] <- "0.1-1%"
	fwrite(dt, file=f, sep="\t", quote=FALSE)
}
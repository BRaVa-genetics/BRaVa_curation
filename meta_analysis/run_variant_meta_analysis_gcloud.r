library(data.table)
library(dplyr)

# Find the METAL files to run, and loop over the meta-analyses
METAL_file_location <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/METAL/variant"
METAL_files <- grep("^.*/METAL_.*.txt$", dir(METAL_file_location, full.names=TRUE), value=TRUE)

METAL_location <- "/well/lindgren/dpalmer/METAL/build/bin/metal"
# Use install_metal.sh script if METAL is not installed

for (file in METAL_files) {
	cat(paste0("carrying out meta-analysis using METAL_file: ", file, "\n"))
	system(paste(METAL_location, file))
}

cat("completed meta-analysis of all traits.\n")

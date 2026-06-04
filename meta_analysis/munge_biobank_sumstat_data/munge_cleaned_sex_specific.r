library(data.table)
library(dplyr)

# Final run through to make sure that sex specific phenotypes (if available) are the only results that are 
# used in the 'cleaned' dataset

# Loop over the folders
# Extract the phenotypes
# Determine if the phenotype is sex-specific
# If it is, and the sex-specific run has been carried out, remove the 'ALL' version

female_sex_specific_phenotypeIDs <- c(
	"Adenomy", "BenCervUterNeo", "BreastCanc", "CervCanc", "Endometr", "FemInf",
	"FemInfAC", "MatHem", "MatHypDis", "OvCanc", "PlacInsuf", "Preeclamps", "UterCanc",
	"PregLoss","PCOS")
male_sex_specific_phenotypeIDs <- c() 

data_folder <- "/Users/dpalmer/Repositories/BRaVa_curation/data/meta_analysis/gcloud"
biobanks <- dir(data_folder)
delete <- TRUE

for (b in biobanks) {
	cat("\n", b, "\n")
	files_gene <- dir(paste0(data_folder, "/", b, "/cleaned/gene"), full.names=TRUE)
	files_variant <- dir(paste0(data_folder, "/", b, "/cleaned/variant"), full.names=TRUE)
	# cat("gene\n")
	# cat(paste0(files_gene, collapse="\n"))
	# cat("\nvariant\n")
	# cat(paste0(files_variant, collapse="\n"))

	for (phe in female_sex_specific_phenotypeIDs)
	{	
		# Gene results files
		if (any(grepl(paste0(phe, "\\.", ".*\\.ALL\\."), files_gene)) & 
			any(grepl(paste0(phe, "\\.", ".*\\.F\\."), files_gene))) {
			cat("Multiple match:\n",
				paste0(grep(paste0(phe, "\\.", ".*\\.ALL\\."), files_gene, value=TRUE),
					collapse='\n'), '\n',
				paste0(grep(paste0(phe, "\\.", ".*\\.F\\."), files_gene, value=TRUE),
					collapse='\n')
				)
			cat("\nRemove the 'ALL' versions, for all ancestries\n")
			for (f in grep(paste0(phe, "\\.", ".*\\.ALL\\."), files_gene, value=TRUE)) {
				if (delete) {
					system(paste("rm", f))
				} else {
					cat(paste("rm", f, "\n"))
				}
			}
			# If there is an ancestry, phenotype match for at least one ancestry,
			# remove all but the sex specific results across the dataset
		} else {
			cat('all good, nothing to remove\n')
		}

		# Variant results files
		if (any(grepl(paste0(phe, "\\.", ".*\\.ALL\\."), files_variant)) & 
			any(grepl(paste0(phe, "\\.", ".*\\.F\\."), files_variant))) {
			cat("Multiple match:\n",
				paste0(grep(paste0(phe, "\\.", ".*\\.ALL\\."), files_variant, value=TRUE),
					collapse='\n'), '\n',
				paste0(grep(paste0(phe, "\\.", ".*\\.F\\."), files_variant, value=TRUE),
					collapse='\n')
				)
			cat("\nRemove the 'ALL' versions, for all ancestries\n")
			for (f in grep(paste0(phe, "\\.", ".*\\.ALL\\."), files_variant, value=TRUE)) {
				if (delete) {
					system(paste("rm", f))
				} else {
					cat(paste("rm", f, "\n"))
				}
			}
			# If there is an ancestry, phenotype match for at least one ancestry,
			# remove all but the sex specific results across the dataset
		} else {
			cat('all good, nothing to remove\n')
		}
	}
}

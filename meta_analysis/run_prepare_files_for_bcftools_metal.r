library(data.table)
library(dplyr)

source("meta_analysis_utils.r")
source("../phenotypes/BRaVa_phenotypes_utils.r")

# Munge the data to run bcftools metal
# export BCFTOOLS_PLUGINS="$HOME/bin"
data_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks"
out_data_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/vcf/variant"

# Place the vcf files
system(paste0("mkdir -p ", out_data_dir))

biobanks <- dir(data_dir)
results_dt_list <- list()
min_cases <- 100

for (biobank in biobanks)
{
	biobank_results_files <- dir(
		paste0(data_dir, "/", biobank, "/cleaned/variant"))
	biobank_results_files_full <- dir(
		paste0(data_dir, "/", biobank, "/cleaned/variant"), full.names=TRUE)
	results <- lapply(biobank_results_files, extract_file_info)
	results_dt_list[[biobank]] <- data.table(
		filename = biobank_results_files_full,
		phenotypeID = sapply(results, `[`, 4),
		pop = sapply(results, `[`, 7),
		binary = sapply(results, `[`, "binary")
		) %>% 
	mutate(n_cases = ifelse(binary, sapply(results, `[`, "n_cases"), NA),
		n_controls = ifelse(binary, sapply(results, `[`, "n_controls"), NA),
		n = ifelse(binary, NA, sapply(results, `[`, "n"))) %>% 
	filter((n > min_cases) | 
		((n_controls > min_cases) & (n_cases > min_cases)))
}

results_dt <- rbindlist(results_dt_list)
phenotypeIDs <- intersect(BRaVa_pilot_phenotypes,
	unique(unlist(results_dt$phenotypeID)))

for (phe in phenotypeIDs) {
	files_variant <- (results_dt %>% filter(phenotypeID == phe))$filename
	for (file in files_variant)  {
		# Read in the variant results file output by SAIGE
		cat(paste(file, "\n"))
		system(paste(
			"sbatch prepare_files_for_bcftools_metal.sh",
			file, out_data_dir))
		break
	}
	break
}

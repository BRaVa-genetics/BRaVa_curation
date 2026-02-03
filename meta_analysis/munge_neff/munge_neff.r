library(data.table)
library(dplyr)
library(ggplot2)

source("../meta_analysis_utils.r")

# Let's combine all of the files together and ensure that the naming conventions are the same across each of the neff files from the biobanks.
# This is all hard-coded, as we have different output formats provided from each biobank
dt_list <- list()

add_ancestry <- function(dt, biobank) {
	if (biobank == "bbj") {
		dt$ancestry <- "EAS"
	} else if (biobank == "genes-and-health") {
		dt$ancestry <- "SAS"
	}
	return(dt)
}

check_pheno <- function(dt) {
	if (any(!(dt$pheno %in% names(file_check_information$phenotype)))) {
		print(dt$pheno[!(dt$pheno %in% names(file_check_information$phenotype))])
		stop(paste0("phenotypes in the Neff file that don't match the expected phenotype names"))
	}
}

check_ancestry <- function(dt) {
	if (any(!(dt$ancestry %in% names(file_check_information$ancestry)))) {
		print(dt$ancestry[!(dt$ancestry %in% names(file_check_information$ancestry))])
		stop(paste0("ancestries in the Neff file that do not match the expected ancestry names"))
	}	
}

renaming_neff_colnames <- list(
	`pheno` = c("pheno", "phenotype", "Phenotype", "Phe", "Pheno"),
	`nglmm` = c("nglmm", "Neff", "nglm"),
	`dataset` = c("dataset", "Dataset", "Biobank", "biobank", "data-set"),
	`ancestry` = c("ancestry", "anc", "pop", "population", "superpop")
)

# Perform checks to ensure that the column names are correct
update_colnames <- function(dt)
{
	for (i in names(renaming_neff_colnames)) {
		print(i)
		if (any(renaming_neff_colnames[[i]] %in% colnames(dt))) {
			colnames(dt)[which(colnames(dt) %in% renaming_neff_colnames[[i]])] <- i
		}
	}
	return(dt)
}

extract_local_neff_indep <- function(folder)
{
	file_info <- list()
	files <- dir(path=folder)
	for (file in files) {
		file_info[[file]] <- extract_file_info(file)
	}
	dt <- rbindlist(file_info, fill=TRUE)
	return(dt)
}

# All of Us
# Here, we need to determine which of these traits to include, and we need to rename them to ensure the 
# phenotype names match those expected by BRaVa.
biobank <- "all-of-us"
neff_filename <- "aou_n_eff_v7.tsv"
dt <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
data_dictionary <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/brava_44_pilot_phenotype_filename.csv"))
names(data_dictionary)[1] <- "phenoname"
names(data_dictionary)[5] <- "ancestry"
setkeyv(data_dictionary, c("phenoname", "ancestry"))
setkeyv(dt, c("phenoname", "ancestry"))
dt <- merge(dt, data_dictionary)
dt <- dt %>% filter(
	((trait_type == "binary") & (phenoname %in% unlist(aou_curated))) |
	(trait_type == "continuous"))
dt <- dt %>% rename(pheno = brava_code, nglmm = n_eff)
# Add in the Neff estimated from the n_cases and n_controls here.
dt$dataset <- "all-of-us"
dt$ancestry <- toupper(dt$ancestry)
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)
dt <- dt %>% mutate(neff_indep = ifelse(is.na(n_controls), n_cases, 4 / ((1/n_cases) + (1/n_controls))))
dt <- dt %>% select(-c("n_cases", "n_controls"))
dt_aou_v7 <- dt
dt_list[[biobank]] <- dt
# Among the same phenotypes, pick the one with the largest number of cases or controls
# (this is what we carried out the first time around).

# Updated All of Us for v8
# This data was generated from the Neff.ipynb in the All of Us researcher workbench
biobank <- "all-of-us"
neff_filename <- "aou_n_eff_v8.tsv.gz"
dt <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename, " | zcat"))
dt <- dt %>% rename(pheno = phenotype, nglmm = Neff)
dt <- dt %>% filter(pheno != "HbA1c")
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)
dt <- dt %>% mutate(neff_indep = ifelse(is.na(n_controls), n_cases, 4 / ((1/n_cases) + (1/n_controls))))
dt <- dt %>% select(-c("n_cases", "n_controls"))
dt_list[[biobank]] <- dt

# Biome
biobank <- "biome"
neff_filename <- "Neff_v2.csv" # Fixed in December 2024
dt <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
dt$dataset <- biobank
dt <- add_ancestry(dt, biobank)
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)

dt_list[[biobank]] <- dt

# BBJ
biobank <- "bbj"
neff_filename <- "pilot_44phenotypes/binary_traits/Nglmm_bbj_binary.20240510.csv"
dt <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
dt$dataset <- biobank
dt <- add_ancestry(dt, biobank)
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)

dt_list[[biobank]] <- dt

# Genes and health
biobank <- "genes-and-health"
neff_filename <- "brava.pilot_phenotypes.Neff.txt"
dt <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
dt$dataset <- biobank
dt <- add_ancestry(dt, biobank)
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)

dt_list[[biobank]] <- dt

# GEL
biobank <- "gel"
# neff_filename <- "pilot_32_phenotypes_20240508/Nglmm/Nglmm_EUR_gel_20240903.csv"
gel_neff_list <- list()
for (anc in c("AFR", "EAS", "EUR", "SAS")) {
	neff_filename <- paste0("pilot_32_phenotypes_allancestries_20250217/Nglmm/Nglmm_", anc, "_v20250217.csv")
	gel_neff_list[[anc]] <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
	gel_neff_list[[anc]]$dataset <- biobank
	# Ancestry is already present
	check_pheno(dt)
	check_ancestry(dt)
}

dt <- rbindlist(gel_neff_list)
dt_list[[biobank]] <- dt

biobank <- "mgbb"
neff_folders_binary <- c("240913_nglmm/nglmm/")
files <- system(paste0("gsutil ls gs://brava-meta-upload-", biobank, "/", neff_folders_binary), intern=TRUE)
# For each of these files, extract the Neff
mgbb_neff_list <- list()
for (file in files) {
	line <- system(paste0("gsutil cat ", file, " | tail -2 | head -1"), intern=TRUE)
	Nglmm <- as.numeric(ifelse(grepl("Nglmm", line), gsub(".* ([0-9\\.]+).*", "\\1", line), NA))
	phenotype <- gsub(".*/([A-Za-z0-9]+)\\..*", "\\1", file)
	anc <- gsub(".*/[A-Za-z0-9]+\\.([A-Z]{3}).*", "\\1", file)
	mgbb_neff_list[[file]] <- data.table(pheno = phenotype, nglmm=Nglmm, dataset="mgbb", ancestry=anc)
	print(mgbb_neff_list[[file]])
}

neff_folders_cts <- c("240913_nglmm/nglmm_qtl/")
# Note that we're hardcoding this separately, because the file convention used for 
# the continuous files differed from the binary files.
# We will also exclude everything but the median value for cts traits.
files <- system(paste0("gsutil ls gs://brava-meta-upload-", biobank, "/", neff_folders_cts), intern=TRUE)
files <- files[which(grepl("_Median\\.", files))]
# For each of these files, extract the Neff
for (file in files) {
	line <- system(paste0("gsutil cat ", file, " | tail -2 | head -1"), intern=TRUE)
	Nglmm <- as.numeric(ifelse(grepl("Nglmm", line), gsub(".* ([0-9\\.]+).*", "\\1", line), NA))
	anc <- gsub(".*/([A-Za-z]{3})\\..*", "\\1", file)
	phenotype <- gsub(".*/[A-Za-z0-9]{3}\\.([A-Za-z0-9]+).*", "\\1", file)
	mgbb_neff_list[[file]] <- data.table(pheno = phenotype, nglmm=Nglmm, dataset="mgbb", ancestry=anc)
}

# Ensure correct renaming of phenotypes
dt <- rbindlist(mgbb_neff_list)
if (any(!(dt$pheno %in% names(file_check_information$phenotype)))) {
	incorrect_name <- dt$pheno[which(!(dt$pheno %in% names(file_check_information$phenotype)))]
	for (name in unique(incorrect_name)) {
		for (n in names(file_check_information$phenotype)) {
			if (name %in% file_check_information$phenotype[[n]]) {
				cat(paste("Renamed:", name, "->", n, "\n"))
				dt$pheno[which(dt$pheno == name)] <- n
				break
			}
		}
	}
}

check_pheno(dt)
check_ancestry(dt)

dt_list[[biobank]] <- dt

# UK Biobank
biobank <- "uk-biobank"
neff_filenames <- c(
	"neff_AFR.csv",
	"neff_AMR.csv",
	"neff_EAS.csv",
	"neff_EUR.csv",
	"neff_SAS.csv"
)

dt_ukb_list <- list()
for (neff_filename in neff_filenames) {
	dt_ukb <- fread(cmd = paste0("gsutil cat gs://brava-meta-upload-", biobank, "/", neff_filename))
	ancestry <- gsub("^neff_([A-Z]{3}).csv", "\\1", neff_filename)
	dt_ukb$ancestry <- ancestry
	dt_ukb_list[[ancestry]] <- dt_ukb
}
dt <- rbindlist(dt_ukb_list)
dt <- dt %>% filter(pheno != "HeightResidSex")
dt$dataset <- biobank
dt <- add_ancestry(dt, biobank)
check_pheno(dt)
check_ancestry(dt)
dt <- update_colnames(dt)

dt_list[[biobank]] <- dt

dt <- rbindlist(dt_list, fill=TRUE)

# Before we write - note that BMI for AMR will be excluded because of inflation, but the v7 version does
# not have that issue. As a result, we can determine Neff from the v7 data, and include that here.
dt <- dt %>% filter(!(dataset == "all-of-us" & ancestry == "AMR" & pheno=="BMI"))
dt <- rbind(dt, dt_aou_v7 %>% filter(ancestry=="AMR", pheno=="BMI"), fill=TRUE)

# fwrite(dt, file = "~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/Neff_weights_may25.tsv.gz", sep='\t', quote=FALSE)
fwrite(dt, file = "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_feb26.tsv.gz", sep='\t', quote=FALSE)
# This file should be placed in the BRaVa_inputs directory: /well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/

dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_may25.tsv.gz")
# Add in column of Neff as estimated by the (Ncases, Ncontrols) or N.
local_gene_sumstat_folders <- c(
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/all-of-us/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/bbj/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/biome/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/genes-and-health/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/gel/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/mgbb/cleaned/gene/",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/uk-biobank/cleaned/gene/"
	)

results_file_information <- list()
for (folder in local_gene_sumstat_folders) {
	results_file_information[[folder]] <- extract_local_neff_indep(folder)
}

results_file_information <- rbindlist(results_file_information, fill=TRUE)
names(results_file_information)[which(names(results_file_information) == "phenotype")] <- "pheno"
setkeyv(results_file_information, c("dataset", "ancestry", "pheno"))
setkeyv(dt, c("dataset", "ancestry", "pheno"))

dt_plot <- merge(dt, results_file_information) %>% mutate(
	n_cases = as.numeric(n_cases),
	n_controls = as.numeric(n_controls),
	n = as.numeric(n)
	)
dt_plot <- dt_plot %>% mutate(neff_ind = ifelse(is.na(n), 4 / ((1/n_cases) + (1/n_controls)), n))

# Merge in the results file information
pdf("../manuscript_figures/Figures/Neff_across_biobanks.pdf", width=10, height=10)
# Plot them against each other, split by (biobank, ancestry)
p <- ggplot(dt_plot) + 
	geom_point(aes(x=neff_ind, y=nglmm, col=ancestry)) + 
	theme_bw() + 
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")  +
	facet_wrap(~ dataset + ancestry, scales="free")
print(p)
dev.off()


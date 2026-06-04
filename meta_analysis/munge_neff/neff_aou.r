library(data.table)
library(dplyr)

results <- list()
i <- 1
for (anc in c("afr", "amr", "eas", "eur", "sas"))
{
	COVAR_LIST <- "age,age2,sex,age_sex,age2_sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20"
	PHENO_FILE <- paste0("allofus_array_", anc ,"_snp_wise_pca_covariates_BMI.csv")
	SPARSE_GRM_FILE <- paste0("allofus_array_", anc, "_snp_wise_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx")
	SPARSE_GRM_ID_FILE <- paste0(SPARSE_GRM_FILE, ".sampleIDs.txt")
	cts_phenotypes <- c(
		"BMI","WHRadjBMI","blood.pressure.diastolic.mean", "blood.pressure.systolic.mean",
		"heart.rate.mean","height","hip.circumference.mean","waist.circumference.mean",
		"weight","median_3004410","median_3006923","median_3007070","median_3013721",
		"median_3020460","median_3022192","median_3027114","median_3028288")

	fwrite(fread(PHENO_FILE) %>% rename(IID = person_id), file="pheno.csv", sep="\t", quote=FALSE)
	PHENO_FILE <- "pheno.csv"

	for (phe in cts_phenotypes) {
		result <- system(paste("Rscript extractNglmm.r",
		        "--phenoFile",  PHENO_FILE,
		        "--phenoCol", phe,
		        "--covarColList", COVAR_LIST,
		        "--traitType", "'quantitative'",
		        "--sparseGRMFile", SPARSE_GRM_FILE,
		        "--sparseGRMSampleIDFile", SPARSE_GRM_ID_FILE,
		        "--useSparseGRMtoFitNULL", TRUE), intern=TRUE
		)
		results[[i]] <- data.table(phenotype = phe, ancestry = anc, Neff = as.numeric(gsub("Nglmm *([0-9.]+) *$", "\\1", result[grep("Nglmm", result)])))
		print(results[[i]])
		i <- i+1
	}

	binary_phenotypes <- c("phe_153","phe_174","phe_174.1","phe_180.1","phe_208","phe_218",
		"phe_250.2","phe_274.1","phe_401","phe_411.2","phe_427.21","phe_440","phe_454",
		"phe_495","phe_531","phe_531.4","phe_550.1","phe_555","phe_585.3","phe_594.1",
		"phe_626.1","phe_626.13","phe_626.14","phe_635","phe_636","phe_636.3","phe_696.4",
		"phe_714.1","phe_747.12","CA_101.41","CV_400.2","CV_401","CV_404.11","CV_416.21",
		"CV_424","CV_431.1","CV_436","CV_440.1","CV_444.1","DE_664.4","EM_202.2","GI_513",
		"GI_520.11","GI_522.1","GI_554.1","GU_585","GU_626.1","GU_626.13","GU_629","MS_703.1",
		"MS_705.1","PP_904.1","RE_474","RE_475","RE_481","SO_374.51","SS_825","phe_4134857",
		"V82")

	for (phe in binary_phenotypes) {
		result <- system(paste("Rscript extractNglmm.r",
		        "--phenoFile",  PHENO_FILE,
		        "--phenoCol", phe,
		        "--covarColList", COVAR_LIST,
		        "--traitType", "'binary'",
		        "--sparseGRMFile", SPARSE_GRM_FILE,
		        "--sparseGRMSampleIDFile", SPARSE_GRM_ID_FILE,
		        "--useSparseGRMtoFitNULL", TRUE), intern=TRUE
		)
		results[[i]] <- data.table(phenotype = phe, ancestry = anc, Neff = as.numeric(gsub("Nglmm *([0-9.]+) *$", "\\1", result[grep("Nglmm", result)])))
		print(results[[i]])
		i <- i+1
	}
}

results <- rbindlist(results)
results$sex <- "ALL"

fwrite(results, file="Neff_ALL.tsv", sep='\t', quote=FALSE)


results <- list()
i <- 1
for (anc in c("afr", "amr", "eas", "eur", "sas"))
{
	COVAR_LIST <- "age,age2,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20"
	PHENO_FILE <- paste0("allofus_array_", anc ,"_snp_wise_pca_covariates_GU_629.csv")
	SPARSE_GRM_FILE <- paste0("allofus_array_", anc, "_snp_wise_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx")
	SPARSE_GRM_ID_FILE <- paste0(SPARSE_GRM_FILE, ".sampleIDs.txt")
	cts_phenotypes <- c(
		"BMI","WHRadjBMI","blood.pressure.diastolic.mean", "blood.pressure.systolic.mean",
		"heart.rate.mean","height","hip.circumference.mean","waist.circumference.mean",
		"weight","median_3004410","median_3006923","median_3007070","median_3013721",
		"median_3020460","median_3022192","median_3027114","median_3028288")

	fwrite(fread(PHENO_FILE) %>% rename(IID = person_id), file="pheno.csv", sep="\t", quote=FALSE)
	PHENO_FILE <- "pheno.csv"

	for (phe in cts_phenotypes) {
		result <- system(paste("Rscript extractNglmm.r",
		        "--phenoFile",  PHENO_FILE,
		        "--phenoCol", phe,
		        "--covarColList", COVAR_LIST,
		        "--traitType", "'quantitative'",
		        "--sparseGRMFile", SPARSE_GRM_FILE,
		        "--sparseGRMSampleIDFile", SPARSE_GRM_ID_FILE,
		        "--useSparseGRMtoFitNULL", TRUE), intern=TRUE
		)
		results[[i]] <- data.table(phenotype = phe, ancestry = anc, Neff = as.numeric(gsub("Nglmm *([0-9.]+) *$", "\\1", result[grep("Nglmm", result)])))
		print(results[[i]])
		i <- i+1
	}

	binary_phenotypes <- c("phe_153","phe_174","phe_174.1","phe_180.1","phe_208","phe_218",
		"phe_250.2","phe_274.1","phe_401","phe_411.2","phe_427.21","phe_440","phe_454",
		"phe_495","phe_531","phe_531.4","phe_550.1","phe_555","phe_585.3","phe_594.1",
		"phe_626.1","phe_626.13","phe_626.14","phe_635","phe_636","phe_636.3","phe_696.4",
		"phe_714.1","phe_747.12","CA_101.41","CV_400.2","CV_401","CV_404.11","CV_416.21",
		"CV_424","CV_431.1","CV_436","CV_440.1","CV_444.1","DE_664.4","EM_202.2","GI_513",
		"GI_520.11","GI_522.1","GI_554.1","GU_585","GU_626.1","GU_626.13","GU_629","MS_703.1",
		"MS_705.1","PP_904.1","RE_474","RE_475","RE_481","SO_374.51","SS_825","phe_4134857",
		"V82")

	for (phe in binary_phenotypes) {
		result <- system(paste("Rscript extractNglmm.r",
		        "--phenoFile",  PHENO_FILE,
		        "--phenoCol", phe,
		        "--covarColList", COVAR_LIST,
		        "--traitType", "'binary'",
		        "--sparseGRMFile", SPARSE_GRM_FILE,
		        "--sparseGRMSampleIDFile", SPARSE_GRM_ID_FILE,
		        "--useSparseGRMtoFitNULL", TRUE), intern=TRUE
		)
		results[[i]] <- data.table(phenotype = phe, ancestry = anc, Neff = as.numeric(gsub("Nglmm *([0-9.]+) *$", "\\1", result[grep("Nglmm", result)])))
		print(results[[i]])
		i <- i+1
	}
}

results <- rbindlist(results)
results$sex <- "F"

fwrite(results, file="Neff_F.tsv", sep='\t', quote=FALSE)

# Extract the collection of traits that were run.
# This can then be passed to the Neff munging script.

# Finally, we can then run meta-analysis over all of the traits (on BMRC).

mapping <- list(
    `BMI` = "BMI",
    `CA_101.41` = "ColonRectCanc",
    `CV_400.2` = "RheumHeaDis",
    `CV_401` = "HTN",
    `CV_404.11` = "CAD",
    `CV_416.21` = "AFib",
    `CV_424` = "HF",
    `CV_431.1` = "Stroke",
    `CV_436` = "PAD",
    `CV_440.1` = "VTE",
    `CV_444.1` = "VaricVeins",
    `DE_664.4` = "Psori",
    `EM_202.2` = "T2Diab",
    `GI_513` = "PeptUlcer",
    `GI_520.11` = "IFHern",
    `GI_522.1` = "IBD",
    `GI_554.1` = "Pancreat",
    `GU_585` = "Urolith",
    `GU_626.1` = "EFRMB",
    `GU_629` = "FemInf",
    `height` = "Height",
    `MS_703.1` = "Gout",
    `MS_705.1` = "RheumArth",
    `PP_904.1` = "MatHem",
    `RE_474` = "COPD",
    `RE_475` = "Asth",
    `RE_481` = "ILDSarc",
    `SO_374.51` = "AMD",
    `WHRadjBMI` = "WHRBMI",
    `phe_174` = "BreastCanc",
    `phe_180.1` = "CervCanc",
    `phe_208` = "BenIntNeo",
    `phe_218` = "BenCervUterNeo",
    `phe_4134857` = "HipRep",
    `phe_585.3` = "CRF",
    `phe_747.12` = "NonRheuValv",
    `median_3007070` = "HDLC",
    `median_3022192` = "TG",
    `median_3027114` = "TChol",
    `median_3028288` = "LDLC",
    `median_3004410` = "HbA1c",
    `median_3006923` = "ALT",
    `median_3013721` = "AST",
    `median_3020460` = "CRP"
)

# Loop through to determine the collection of outputs to merge on
files <- system("ls ~/conditioning_test/summary_stats_ensure_no_PID", intern=TRUE)
# Remove the files that are the subset of EFRMB with fewer cases
files <- setdiff(files,
	c("all-of-us.palmer.BRaVa_v8.EFRMB.v8.F.AFR.1898.32309.SAIGE.gene.20240326.txt.gz",
	  "all-of-us.palmer.BRaVa_v8.EFRMB.v8.F.AMR.1785.31143.SAIGE.gene.20240326.txt.gz",
	  "all-of-us.palmer.BRaVa_v8.EFRMB.v8.F.EAS.197.3711.SAIGE.gene.20240326.txt.gz",
	  "all-of-us.palmer.BRaVa_v8.EFRMB.v8.F.EUR.12234.84572.SAIGE.gene.20240326.txt.gz")
)

extract_file_info <- function(filename)
{
	gz <- ifelse(grepl(".gz$", filename), TRUE, FALSE)
	filename <- gsub(".gz$", "", filename)
	filename <- gsub("cleaned.", "", filename)
	file_info <- as.list(strsplit(filename, split="\\.")[[1]])

	if (!(length(file_info) %in% c(12, 13))) {
		print(file_info)
		stop("Incorrect file naming convention, please check filenames")
	}

	if (length(file_info) == 13) {
		binary <- TRUE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n_cases",
			"n_controls",
			"software",
			"type",
			"date",
			"split")
	} else {
		binary <- FALSE
		names(file_info) <- c(
			"dataset",
			"last_name",
			"analysis_name",
			"phenotype",
			"freeze_number",
			"sex",
			"ancestry",
			"n",
			"software",
			"type",
			"date",
			"split")
	}

	file_info$gz <- gz
	file_info$binary <- binary
	return(file_info)
}

dt_files <- list()
for (file in files) {
	dt_files[[file]] <- extract_file_info(file)
}
dt_files <- rbindlist(dt_files, fill=TRUE)
setkeyv(dt_files, c("phenotype", "ancestry", "sex"))

dt <- fread("Neff_ALL.tsv")
dt <- dt %>% filter(phenotype %in% names(mapping))
dt <- dt %>% mutate(phenotype = unlist(mapping[dt$phenotype]), ancestry=toupper(ancestry))
setkeyv(dt, c("phenotype", "ancestry", "sex"))

dt_all <- merge(dt_files, dt)

dt <- fread("Neff_F.tsv")
dt <- dt %>% filter(phenotype %in% names(mapping))
dt <- dt %>% mutate(phenotype = unlist(mapping[dt$phenotype]), ancestry=toupper(ancestry))
setkeyv(dt, c("phenotype", "ancestry", "sex"))

dt_f <- merge(dt_files, dt)

dt <- rbind(dt_all, dt_f)

fwrite(dt, file="aou_neff_v8.tsv.gz", sep="\t", quote=FALSE)




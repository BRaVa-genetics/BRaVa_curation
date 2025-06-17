library(data.table)
library(dplyr)
library(ggplot2)
library(googlesheets4)
library(stringr)

BRaVa_pilot_phenotypes <- c(
	"AMD", "Asth", "AFib", "BenCervUterNeo", "BenIntNeo",
	"BreastCanc", "CervCanc", "COPD", "CRF", "ColonRectCanc",
	"CAD", "EFRMB", "FemInf", "Gout", "HF", "HTN", "IBD",
	"IFHern",    "ILDSarc",   "MatHem", "NonRheuValv", "Pancreat",
	"PeptUlcer", "PAD", "Psori", "RheumHeaDis", "RheumArth",
	"Stroke", "T2Diab", "Urolith", "VaricVeins", "VTE", "ALT",
	"AlcCons", "AST", "BMI", "CRP", "HDLC", "Height", "LDLC", 
	"TChol", "TG", "WHRBMI", "HipRep")

extract_BRaVa_pilot_phenotypes <- function(pilot_only=TRUE) {
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249",
		sheet="Sequenced_Sample_Sizes", skip=5)
	dt_pilot <- dt %>% select("Phenotype ID", "(ICD_Phecode/Procedures tab)", "...80", "...113")
	names(dt_pilot) <- c("phenotypeID", "case_control", "count", "include")
	dt_pilot <- dt_pilot %>% mutate(include = as.logical(include))
	if (pilot_only) { dt_pilot <- dt_pilot %>% filter(include) }
	dt_pilot <- dt_pilot %>% filter(phenotypeID != "")
	return(dt_pilot$phenotypeID)
}

extract_BRaVa_pilot_phenotypes_ukb_counts <- function(pilot_only=TRUE) {
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249",
		sheet="Sequenced_Sample_Sizes", skip=5)
	dt_pilot <- dt %>% select("Phenotype ID", "(ICD_Phecode/Procedures tab)", "...80", "...113")
	names(dt_pilot) <- c("phenotypeID", "case_control", "count", "include")
	dt_pilot <- dt_pilot %>% mutate(include = as.logical(include))
	if (pilot_only) { dt_pilot <- dt_pilot %>% filter(include) }
	dt_pilot <- dt_pilot %>% filter(phenotypeID != "")
	return(dt_pilot %>% select(phenotypeID, case_control, count))
}

munge_BRaVa_ICD_proposals <- function()
{
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet="ICD_Phecode", skip=3)
	cols <- c(
		"Description", "Phenotype ID",
		"ICD10_control_exclude", "ICD10_case_include",
		"ICD9_control_exclude", "ICD9_case_include")
	dt <- dt[, cols, with=FALSE]
	
	# Remove empty rows
	dt <- dt %>% filter(!is.na(Description))
	names(dt) <- c("phenotype", "phenotypeID", "ICD10_control_exclude", "ICD10_case_include", "ICD9_control_exclude", "ICD9_case_include")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, ICD10_control_exclude := gsub("[[:space:]]", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("[[:space:]]", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("[[:space:]]", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("[[:space:]]", "", ICD9_case_include)]
	
	# Replace '.'s with '' (it's an equivalent encoding for ICD codes).
	dt[, ICD10_control_exclude := gsub("\\.", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("\\.", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("\\.", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("\\.", "", ICD9_case_include)]

	# Replace commas for pipes to get the regular expression ready.
	dt[, ICD10_control_exclude := gsub(",|;", " | ", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub(",|;", " | ", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub(",|;", " | ", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub(",|;", " | ", ICD9_case_include)]

	# Ensure that the first and final character of the regular expression is a space.
	dt[, ICD10_control_exclude := gsub("^(.*)$", " \\1 ", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("^(.*)$", " \\1 ", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("^(.*)$", " \\1 ", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("^(.*)$", " \\1 ", ICD9_case_include)]

	# Remove asterisks
	dt[, ICD10_control_exclude := gsub("\\*", "", ICD10_control_exclude)]
	dt[, ICD10_case_include := gsub("\\*", "", ICD10_case_include)]
	dt[, ICD9_control_exclude := gsub("\\*", "", ICD9_control_exclude)]
	dt[, ICD9_case_include := gsub("\\*", "", ICD9_case_include)]

	dt$ICD10_control_exclude[which(is.na(dt$ICD10_control_exclude))] <- ""
	dt$ICD10_case_include[which(is.na(dt$ICD10_case_include))] <- ""
	dt$ICD9_control_exclude[which(is.na(dt$ICD9_control_exclude))] <- ""
	dt$ICD9_case_include[which(is.na(dt$ICD9_case_include))] <- ""

	return(dt)
}

munge_BRaVa_cts_ICD_exclusion_proposals <- function()
{
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249",
		sheet="ICD_Phecode_cts_exclusion")
	cols <- c("Description", "Phenotype ID", "ICD10_exclude", "ICD9_exclude")
	dt <- dt[, cols, with=FALSE]

	# Remove empty rows
	dt <- dt %>% filter(!is.na(Description))
	names(dt) <- c("phenotype", "phenotypeID", "ICD10_exclude", "ICD9_exclude")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, ICD10_exclude := gsub("[[:space:]]", "", ICD10_exclude)]
	dt[, ICD9_exclude := gsub("[[:space:]]", "", ICD9_exclude)]
	
	# Replace '.'s with '' (it's an equivalent encoding for ICD codes).
	dt[, ICD10_exclude := gsub("\\.", "", ICD10_exclude)]
	dt[, ICD9_exclude := gsub("\\.", "", ICD9_exclude)]

	# Replace commas for pipes to get the regular expression ready.
	dt[, ICD10_exclude := gsub(",|;", " | ", ICD10_exclude)]
	dt[, ICD9_exclude := gsub(",|;", " | ", ICD9_exclude)]

	# Ensure that the first and final character of the regular expression is a space.
	dt[, ICD10_exclude := gsub("^(.*)$", " \\1 ", ICD10_exclude)]
	dt[, ICD9_exclude := gsub("^(.*)$", " \\1 ", ICD9_exclude)]

	# Remove asterisks
	dt[, ICD10_exclude := gsub("\\*", "", ICD10_exclude)]
	dt[, ICD9_exclude := gsub("\\*", "", ICD9_exclude)]

	dt$ICD10_exclude[which(is.na(dt$ICD10_exclude))] <- ""
	dt$ICD9_exclude[which(is.na(dt$ICD9_exclude))] <- ""

	return(dt)
}

extract_ID_and_ICD_UKB <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv"
)
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Extract the relevant columns for each of the encodings
	ICD10s <- c("41202", "41204", "40006", "40001", "40002")
	ICD9s <- c("41203", "41205", "40013")

	cols <-  get_cols(c(ICD9s, ICD10s), dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)
	
	# Extract the ICD10 columns
	ICD10_cols <-  get_cols(ICD10s, dt)

	# Combine vector into a single string
	dt[, ICD10_string := do.call(paste,.SD), .SDcols=ICD10_cols]
	# Remove NAs
	dt[, ICD10_string := gsub("NA", "", ICD10_string)]
	# Pad start and end with space
	dt[, ICD10_string := gsub("^(.*)$", " \\1 ", ICD10_string)]
	# Remove trailing whitespace
	dt[, ICD10_string := gsub("( )+", " ", ICD10_string)]
	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, ICD10_string := gsub("\\.", "", ICD10_string)]

	# Extract the ICD9 columns
	ICD9_cols <-  get_cols(ICD9s, dt)

	# Combine vector into a single string
	dt[, ICD9_string := do.call(paste,.SD), .SDcols=ICD9_cols]
	# Remove NAs
	dt[, ICD9_string := gsub("NA", "", ICD9_string)]
	# Pad start and end with space
	dt[, ICD9_string := gsub("^(.*)$", " \\1 ", ICD9_string)]
	# Remove trailing whitespace
	dt[, ICD9_string := gsub("( )+", " ", ICD9_string)]

	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, ICD9_string := gsub("\\.", "", ICD9_string)]

	if (length(which(is.na(dt$ICD10_string))) > 0) {
		dt$ICD10_string[which(is.na(dt$ICD10_string))] <- ""
	}

	if (length(which(is.na(dt$ICD9_string))) > 0) {
		dt$ICD9_string[which(is.na(dt$ICD9_string))] <- ""
	}

	cols <- c("eid", "ICD9_string", "ICD10_string")
	dt <- dt[, cols, with=FALSE]
	return(dt)
}

simplify_query <- function(query)
{
	# Simplify to speed things up
	query_tmp <- gsub("^ (.*) ", "\\1", query)
	query_tmp <- strsplit(query_tmp, split=" \\| ")[[1]]
	query_tmp <- query_tmp[order(nchar(query_tmp))]
	query <- c()
	while (length(query_tmp) > 0) {
		query <- c(query, query_tmp[1])
		match <- grep(query[length(query)], query_tmp)
		query_tmp <- query_tmp[-match]

	}
	query <- paste0(query, collapse="|(^|\\s)")
	query <- paste0("(^|\\s)", query)
	return(query)
}

extract_regexp <- function(dt, cols)
{
	for (col in cols) {
		for (i in 1:nrow(dt)) {
			current_regexp <- dt[[col]][i]
			dt[[col]][i] <- ifelse(
				current_regexp == "", "", simplify_query(current_regexp))
		}
	}
	return(dt)
}

extract_case_status <- function(dt_data, dt_query,
	assume_tree=TRUE, just_tree_query=FALSE) {
	
	cols <- dt_query %>% select(phenotype, phenotypeID)
	
	if (!just_tree_query) {
		dt_data[, (cols$phenotypeID) := 0]
	}

	for (i in 1:nrow(cols)) {
		# Add a new column to dt_data
		phenotypeID <- cols$phenotypeID[i]
		cat(paste0(cols$phenotype[i], " (phenotype ID: ", phenotypeID, ")...\n"))

		if (dt_query$ICD10_control_exclude[i] != "")
		{
			query <- dt_query$ICD10_control_exclude[i]
			if (assume_tree) {
				query <- simplify_query(query)
				cat(paste0("query regular expression: ", query, "\n"))
			}

			if (!just_tree_query) {
				ICD10_where <- grep(query, dt_data$ICD10_string)
				cat(paste0("ICD10 control exclusions: ", length(ICD10_where), "...\n"))
				dt_data[[phenotypeID]][ICD10_where] <- NA
			}
			dt_query$ICD10_control_exclude[i] <- query
		}

		if (dt_query$ICD9_control_exclude[i] != "")
		{
			query <- dt_query$ICD9_control_exclude[i]
			if (assume_tree) {
				query <- simplify_query(query)
				cat(paste0("query regular expression: ", query, "\n"))
			}

			if (!just_tree_query) {
				ICD9_where <- grep(query, dt_data$ICD9_string)
				cat(paste0("ICD9 control exclusions: ", length(ICD9_where), "...\n"))
				dt_data[[phenotypeID]][ICD9_where] <- NA
			}
			dt_query$ICD9_control_exclude[i] <- query
		}

		if (!just_tree_query) {
			cat(paste0(
				"total control exclusions count prior to adding cases: ",
				sum(is.na(dt_data[[phenotypeID]])), "\n"))
		}

		if (dt_query$ICD10_case_include[i] != "")
		{
			query <- dt_query$ICD10_case_include[i]
			if (assume_tree) {
				query <- simplify_query(query)
				cat(paste0("query regular expression: ", query, "\n"))
			}

			if (!just_tree_query) {
				ICD10_where <- grep(query, dt_data$ICD10_string)
				cat(paste0("ICD10 cases: ", length(ICD10_where), "...\n"))
				dt_data[[phenotypeID]][ICD10_where] <- 1
			}
			dt_query$ICD10_case_include[i] <- query
		}

		if (dt_query$ICD9_case_include[i] != "")
		{
			query <- dt_query$ICD9_case_include[i]
			if (assume_tree) {
				query <- simplify_query(query)
				cat(paste0("query regular expression: ", query, "\n"))
			}

			if (!just_tree_query) {
				ICD9_where <- grep(query, dt_data$ICD9_string)
				cat(paste0("ICD9 cases: ", length(ICD9_where), "...\n"))
				dt_data[[phenotypeID]][ICD9_where] <- 1
			}
			dt_query$ICD9_case_include[i] <- query
		}

		if (!just_tree_query) {
			cat(paste0("total case count: ", sum(dt_data[[phenotypeID]], na.rm=TRUE), "\n"))
			cat(paste0("final case count: ", sum(dt_data[[phenotypeID]], na.rm=TRUE), "\n"))
			cat(paste0("final control exclusions count after adding cases: ",
				sum(is.na(dt_data[[phenotypeID]])), "\n"))
			cat(paste0("final control count: ",
				sum(dt_data[[phenotypeID]] == 0, na.rm=TRUE), "\n\n"))
		}
	}

	return(list(dt_data=dt_data, dt_query=dt_query))
}

extract_continuous_trait_counts <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv",
	biomarker_file = "/well/lindgren/UKBIOBANK/DATA/Biomarker_data/ukb27722.csv",
	manual_curation_mapping = "data/BRaVa_continuous_trait_manual_mappings.tsv",
	superpopulation_labels = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/superpopulation_labels.tsv",
	write_continuous_data_file = TRUE,
	counts_out = "data/output",
	file_out = "/well/lindgren/UKBIOBANK/dpalmer/superpopulation_assignments/BRaVa_cts_phenotypes_with_superpopulation_labels_updated.tsv.gz"
)
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { 
			cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE))
		}
		return(cols)
	}

	dt_pheno <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# First the columns in the phenotype file
	dt_manual <- fread(manual_curation_mapping) %>% mutate(UKB_code = as.character(UKB_code))
	cols <- (dt_manual %>% filter(UKB_code != ""))$UKB_code
	cols <-  get_cols(cols, dt_pheno)
	pheno_cols <- grep("*\\-0\\.0", cols, value=TRUE)
	select_cols <- c("character", rep("numeric", length(pheno_cols)))
	names(select_cols) <- c("eid", pheno_cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt_pheno <- fread(phenotype_file, select=select_cols)
	setkey(dt_pheno, "eid")

	# Next, the columns in the biomarker file
	dt_biomarker <- fread(biomarker_file, na.strings=NULL, nrow=1)

	cols <- (dt_manual %>% filter(UKB_code != ""))$UKB_code
	cols <-  get_cols(cols, dt_biomarker)
	biomarker_cols <- grep("*\\-0\\.0", cols, value=TRUE)
	select_cols <- c("character", rep("numeric", length(biomarker_cols)))
	names(select_cols) <- c("eid", biomarker_cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt_biomarker <- fread(biomarker_file, select=select_cols)

	# Get population specific counts
	# Merge with 1000G labels - this file is created using 05_estimate_superpopulation.r in the QC folder.
	dt_classify <- fread(superpopulation_labels)
	names(dt_classify) <- gsub("^PC([0-9]+)$", "kg_pc\\1", names(dt_classify))
	dt_classify[, eid:=sample.ID]
	dt_classify[, sample.ID:=NULL]
	setkey(dt_classify, "eid")

	dt_cts_classified <- merge(merge(dt_pheno, dt_classify), dt_biomarker)
	
	dt_cts_classified <- dt_cts_classified %>% mutate(WHR = (`48-0.0` / `49-0.0`))
	model <- lm(WHR ~ `21001-0.0`, data=dt_cts_classified)
	WHRBMI <- data.table(
		eid = dt_cts_classified$eid[as.integer(names(resid(model)))], WHRBMI = resid(model))
	setkey(WHRBMI, "eid")
	dt_cts_classified <- merge(dt_cts_classified, WHRBMI, all=TRUE)
	
	# Map the names back
	tmp <- names(dt_cts_classified)
	tmp <- gsub("-.*", "", tmp)
	new_phenocol <- c()
	for (i in 1:length(tmp)) {
		if (tmp[i] %in% dt_manual$UKB_code) {
			tmp[i] <- dt_manual$phenotypeID[which(dt_manual$UKB_code == tmp[i])]
			new_phenocol <- c(new_phenocol, tmp[i])
		}
	}
	names(dt_cts_classified) <- tmp

	if (write_continuous_data_file) {
		fwrite(dt_cts_classified, file=file_out, sep="\t")
	}

	sum_not_is.na <- function(col) { sum(!is.na(col)) }

	# Split by 1000G label and count the number of non-NA entries
	dt_counts <- dt_cts_classified %>% 
		group_by(classification_strict) %>% 
		summarise(across(c(new_phenocol, "WHRBMI"), sum_not_is.na))
	dt_counts_t <- data.table(data.table::transpose(
		dt_counts, keep.names="phenotypeID", make.names="classification_strict"))
	setkey(dt_counts_t, "phenotypeID")
	setkey(dt_manual, "phenotypeID")
	dt_counts_t <- merge(dt_manual, dt_counts_t, all.y=TRUE)
	dt_counts_t <- dt_counts_t %>% select(-phenotype)
	fwrite(dt_counts_t, file=paste0(counts_out, "/UKBB_cts_non_missing_counts.tsv.gz"), sep="\t")

	dt_counts <- dt_cts_classified %>% summarise(across(c(new_phenocol, "WHRBMI"), sum_not_is.na))
	dt_counts_t <- data.table::transpose(dt_counts, keep.names="phenotypeID")
	dt_counts_t <- merge(dt_manual, dt_counts_t, all.y=TRUE)
	dt_counts_t <- dt_counts_t %>% select(-phenotype)

	fwrite(dt_counts_t, file=paste0(counts_out, "/UKBB_cts_non_missing_counts_total.tsv.gz"), sep="\t")
}

munge_BRaVa_OPCS4_proposals <- function() {
	# Download and extract the case and control ICD9/10 codes from the BRaVa nominate phenotypes file
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249", sheet="Procedures")
	cols <- c("Procedures", "Phenotype ID", "OPCS codes")

	dt <- dt[, cols, with=FALSE]
	names(dt) <- c("phenotype", "phenotypeID", "OPCS4_code")
	dt <- data.table(dt)

	# Remove any extra spaces
	dt[, OPCS4_code := gsub("[[:space:]]", "", OPCS4_code)]
	# Replace '.'s with '' (it's an equivalent encoding for OPCS4 codes).
	dt[, OPCS4_code := gsub("\\.", "", OPCS4_code)]
	# Replace commas and semi-colons for pipes to get the regular expression ready.
	dt[, OPCS4_code := paste0(
		"(^|\\s)",
		gsub(",|;", "(\\\\s|$)|(^|\\\\s)", OPCS4_code),
		"(\\s|$)")]
	# Remove asterisks
	dt[, OPCS4_code := gsub("\\*", "", OPCS4_code)]
	dt$OPCS4_code[which(is.na(dt$IOPCS4_code))] <- ""

	return(dt)
}

extract_ID_and_OPCS4_UKB <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv")
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Extract the relevant columns for each of the encodings
	OPCS4s <- c("41200", "41210", "41272")
	cols <-  get_cols(OPCS4s, dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)

	# Combine vector into a single string
	dt[, OPCS4_string := do.call(paste,.SD), .SDcols=cols]
	# Remove NAs
	dt[, OPCS4_string := gsub("NA", "", OPCS4_string)]
	# Pad start and end with space
	dt[, OPCS4_string := gsub("^(.*)$", " \\1 ", OPCS4_string)]
	# Remove trailing whitespace
	dt[, OPCS4_string := gsub("( )+", " ", OPCS4_string)]
	# Remove the '.'s from the file (there aren't any, but this will be useful for collaborators).
	dt[, OPCS4_string := gsub("\\.", "", OPCS4_string)]

	if (length(which(is.na(dt$OPCS4_string))) > 0) {
		dt$OPCS4_string[which(is.na(dt$OPCS4_string))] <- ""
	}

	cols <- c("eid", "OPCS4_string")
	dt <- dt[, cols, with=FALSE]
	return(dt)
}

extract_procedure_case_status <- function(dt_data, dt_query)
{	
	cols <- dt_query %>% select(phenotype, phenotypeID)
	dt_data[, (cols$phenotypeID) := 0]
	
	for (i in 1:nrow(cols)) {
		# Add a new column to dt_data
		phenotypeID <- cols$phenotypeID[i]
		cat(paste0(cols$phenotype[i], " (phenotype ID: ", phenotypeID, ")...\n"))

		if (dt_query$OPCS4_code[i] != "")
		{
			query <- dt_query$OPCS4_code[i]
			cat(paste0("query regular expression: ", query, "\n"))
			OPCS4_where <- grep(query, dt_data$OPCS4_string)
			cat(paste0("OPCS4 code cases: ", length(OPCS4_where), "...\n"))
			dt_data[[phenotypeID]][OPCS4_where] <- 1
		}
		cat(paste0("total case count: ", sum(dt_data[[phenotypeID]], na.rm=TRUE), "\n"))
	}
	return(dt_data)
}

extract_covariates <- function(
	phenotype_file = "/well/lindgren-ukbb/projects/ukbb-11867/DATA/PHENOTYPE/PHENOTYPE_MAIN/ukb10844_ukb50009_updateddiagnoses_14012022.csv")
{
	get_cols <- function(codes, dt, na.filter=FALSE)
	{
		cols <- c()
		for (code in codes) { cols <- c(cols, grep(paste0("^", code, "\\-"), names(dt), value=TRUE)) }
		return(cols)
	}

	dt <- fread(phenotype_file, na.strings=NULL, nrow=1)

	# Extract the relevant columns for each of the encodings
	covs <- c("22001", "21022", "22009")
	cols <-  get_cols(covs, dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)
	names(dt) <- c("eid", "sex", "age", paste0("ukbb_pc", seq(1, sum(grepl("^22009-", names(dt))))))
	return(dt)
}

extract_terms <- function(
    strings_to_map,
    mapping_file = "~/Repositories/BRaVa_curation/phenotypes/data/all_lkps_maps_v4.xlsx",
    mapping_file_sheet = "icd9_icd10",
    from_col="ICD10",
    to_col="ICD9",
    grep=FALSE)
{
    # Read in the excel spreadsheet converting x -> y
    # This includes all possible x terms and what they're mapped to.
    dt <- data.table(read_excel(mapping_file, sheet = mapping_file_sheet))
    if (grep) {
    	matches <- grep(strings_to_map, dt[[from_col]])
    	if (length(matches) > 0) {
    		return(setdiff(unique(dt[[to_col]][matches]), "UNDEF"))
		} else {
			return(c())
		}
	} else {
		return(setdiff(unique(dt[get(from_col) %in% strings_to_map][,get(to_col)]), "UNDEF"))
	}
}



# The diseases and their prevalences.
library(xml2)
library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)

# Permalink for the query
# https://ghdx.healthdata.org/gbd-results-tool?params=gbd-api-2019-permalink/4a29a063ee2f584a730e4df000fbd9e7

dt <- fread("data/IHME-GBD_2019_DATA-1b2174ce-1.csv")
dt <- dt %>% filter(age_name != "All Ages", sex_name == "Both") %>% filter(metric_name == "Percent", val >= 0.00025, val <= 0.01) %>% group_by(cause_name) %>% summarise(max_prev = max(val))

pdf(file="plot_rareish_GBD.pdf", width=10, height=10)
p <- ggplot(dt %>% mutate(cause_name=reorder(cause_name, max_prev)), aes(x=max_prev*100000, y=cause_name)) + 
	geom_point() + 
	scale_x_continuous(trans='log10') +
	xlab("Point prevalence (per 100,000)") + 
	ylab("") + 
	theme_classic() + theme(axis.text.y = element_text(size = 6)) +
	geom_vline(xintercept = 25, color="red") + 
	geom_vline(xintercept = 1000, color="red")
print(p)
dev.off()

fwrite(dt, file="data/output/GBD_summarised_adult_diseases.tsv", sep='\t')
# This file is then pasted into google sheets, and manually curated (sheet 1)

ukb_icd_counts <- function(
	dt_ICD_codes_to_lookup,
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
	# Remove trailing whitespace
	dt[, ICD10_string := gsub("( )+", " ", ICD10_string)]
	# Pad start and end with space
	dt[, ICD10_string := gsub("^(.*)$", " \\1 ", ICD10_string)]

	# Extract the ICD9 columns
	ICD9_cols <-  get_cols(ICD9s, dt)

	# Combine vector into a single string
	dt[, ICD9_string := do.call(paste,.SD), .SDcols=ICD9_cols]
	# Remove trailing whitespace
	dt[, ICD9_string := gsub("( )+", " ", ICD9_string)]
	# Pad start and end with space
	dt[, ICD9_string := gsub("^(.*)$", " \\1 ", ICD9_string)]

	dt_ICD_codes_to_lookup$count <- NA
	for (i in 1:nrow(dt_ICD_codes_to_lookup)) {
		print(dt_ICD_codes_to_lookup$Cause[i])
		grep_exp_ICD10 <- dt_ICD_codes_to_lookup$coding_ICD10_string[i]
		grep_exp_ICD9 <- dt_ICD_codes_to_lookup$coding_ICD9_string[i]
		if (grep_exp_ICD10 == "") {
			if (grep_exp_ICD9 == "") {
				print("Both empty")
				dt_ICD_codes_to_lookup$count[i] <- 0
			} else {
				print("ICD10 empty, ICD9 not empty")
				dt_ICD_codes_to_lookup$count[i] <- sum(grepl(grep_exp_ICD9, dt$ICD9_string))
			}
		} else {
			# ICD10 is not empty
			if (grep_exp_ICD9 == "") {
				print("ICD9 empty, ICD10 not empty")
				dt_ICD_codes_to_lookup$count[i] <- sum(grepl(grep_exp_ICD10, dt$ICD10_string))
			} else {
				print("ICD9 and ICD10 not empty")
				dt_ICD_codes_to_lookup$count[i] <- sum((grepl(grep_exp_ICD10, dt$ICD10_string) | grepl(grep_exp_ICD9, dt$ICD9_string)))
			}
		}
		print(dt_ICD_codes_to_lookup$count[i])
	}
	dt_ICD_codes_to_lookup <- data.table(dt_ICD_codes_to_lookup %>% select(Cause, coding_ICD10_string, coding_ICD9_string, count))
	setkey(dt_ICD_codes_to_lookup, "Cause")
	return(dt_ICD_codes_to_lookup)
}

create_GBD_list <- function(
	GBD_icd_codes="data/icd_mappings/IHME_GBD_2019_CAUSE_ICD_CODE_MAP_Y2020M10D15.XLSX",
	skip_GBD=1) {
	dt <- read_excel(GBD_icd_codes, skip=skip_GBD) %>% select(Cause, ICD10, ICD9)
	dt <- dt %>% mutate(
		ICD10=gsub("[[:space:]]", "", ICD10),
		ICD9=gsub("[[:space:]]", "", ICD9),
		cause_name=Cause
		) %>% mutate(ICD10=gsub("\\.", "", ICD10), ICD9=gsub("\\.", "", ICD9))
	return(data.table(dt, key="cause_name"))
}

dt <- fread("data/output/GBD_summarised_adult_diseases.tsv", key="cause_name")
dt <- merge(dt, create_GBD_list())
# Now extract all of these ICD codes

# Extract all the ICD codes that are single ICD codes
replace_NA_empty <- function(vector) { 
	vector[which(is.na(vector))] <- "" 
	return(vector)
}

dt_ICD10_list <- lapply(strsplit(dt$ICD10, split=","), grep, pattern="\\-", value=TRUE, invert=TRUE)
dt_ICD10_list <- lapply(dt_ICD10_list, replace_NA_empty)
dt$ICD10_string <- unlist(lapply(dt_ICD10_list, paste, collapse=","))

dt_ICD10_list <- lapply(strsplit(dt$ICD10, split=","), grep, pattern="\\-", value=TRUE)
dt_ICD10CM_ordered <- gsub("\\.", "", gsub("^.*<name>(.*)</name>", "\\1", grep(".*<name>.*</name>", readLines("data/icd_mappings/cms_gov_files/ICD10/icd10cm_tabular_2021.xml"), value=TRUE)))
dt_ICD10CM <- grep("^[0-9]", dt_ICD10CM_ordered, invert=TRUE, value=TRUE)
dt_ICD10CM <- data.table(coding = sort(dt_ICD10CM))
setkey(dt_ICD10CM, key="coding")
dt_ICD10_UKB <- fread("data/icd_mappings/coding19.tsv", key="coding") # These are all in order. The filtering suggested in the spreadsheet assumes alphanumeric order.
dt_ICD10 <- merge(dt_ICD10CM, dt_ICD10_UKB, all=TRUE) %>% select(coding)

# Read in the ICD10 codes
i <- 1
for (ICD_codes in dt_ICD10_list) {
	for (ICD_code in ICD_codes) {
		ICDs <- unlist(strsplit(ICD_code, split="-"))
		start <- which(dt_ICD10$coding == ICDs[1])
		stop <- which(dt_ICD10$coding == ICDs[2])
		if ((length(start) == 0) | (length(stop) == 0)) {
			if (length(stop) == 0) {
				print(paste("stop is empty", i, ICDs))
				if (length(start) == 0) {
					print("damn")
					next
				} else {
					j <- 0
					while(sort(c(dt_ICD10$coding[start+j], ICDs[2]))[2] == ICDs[2]) {
						j <- j+1
					}
					stop <- start+j-1
				}
			}
			if (length(start) == 0) {
				print(paste("start is empty", i, ICDs))
				j <- 0
				while(sort(c(dt_ICD10$coding[stop+j], ICDs[1]))[1] == ICDs[1]) {
					j <- j-1
				}
				start <- stop+j+1
			}
		}
		ICDs_to_add <- dt_ICD10$coding[seq(start, stop)]
		ICDs_to_add <- intersect(ICDs_to_add, dt_ICD10_UKB$coding)
		if (dt$ICD10_string[i] != "") {
			dt$ICD10_string[i] <- paste(dt$ICD10_string[i], paste(ICDs_to_add, collapse=","), sep=",")
		} else {
			dt$ICD10_string[i] <- paste(ICDs_to_add, collapse=",")
		}
	}
	i <- i+1
}

dt$ICD10_string <- gsub("\\,+", "\\,", dt$ICD10_string)
dt$ICD10_string <- gsub("^\\,", "", dt$ICD10_string)
dt$ICD10_string <- gsub("\\,$", "", dt$ICD10_string)

dt_ICD9_list <- lapply(strsplit(dt$ICD9, split=","), grep, pattern="\\-", value=TRUE, invert=TRUE)
dt_ICD9_list <- lapply(dt_ICD9_list, replace_NA_empty)
dt$ICD9_string <- unlist(lapply(dt_ICD9_list, paste, collapse=","))

dt_ICD9_list <- lapply(strsplit(dt$ICD9, split=","), grep, pattern="\\-", value=TRUE)
dt_ICD9CM <- gsub(" .*", "", readLines("data/icd_mappings/cms_gov_files/ICD9/CMS32_DESC_LONG_DX.txt"))
dt_ICD9CM <- data.table(coding = sort(dt_ICD9CM))
setkey(dt_ICD9CM, key="coding")
dt_ICD9_UKB <- fread("data/icd_mappings/coding87.tsv", key="coding") # These are all in order. The filtering suggested in the spreadsheet assumes alphanumeric order.
dt_ICD9 <- merge(dt_ICD9CM, dt_ICD9_UKB, all=TRUE)

# Read in the ICD9 codes
i <- 1
for (ICD_codes in dt_ICD9_list) {
	for (ICD_code in ICD_codes) {
		ICDs <- unlist(strsplit(ICD_code, split="-"))
		start <- which(dt_ICD9$coding == ICDs[1])
		stop <- which(dt_ICD9$coding == ICDs[2])
		if ((length(start) == 0) | (length(stop) == 0)) {
			if (length(stop) == 0) {
				print(paste("stop is empty", i, ICDs))
				if (length(start) == 0) {
					print("damn")
					next
				} else {
					j <- 0
					while(sort(c(dt_ICD9$coding[start+j], ICDs[2]))[2] == ICDs[2]) {
						j <- j+1
					}
					stop <- start+j-1
				}
			}
			if (length(start) == 0) {
				print(paste("start is empty", i, ICDs))
				j <- 0
				while(sort(c(dt_ICD9$coding[stop+j], ICDs[1]))[1] == ICDs[1]) {
					j <- j-1
				}
				start <- stop+j+1
			}
		}
		ICDs_to_add <- dt_ICD9$coding[seq(start, stop)]
		ICDs_to_add <- intersect(ICDs_to_add, dt_ICD9_UKB$coding)
		if (dt$ICD9_string[i] != "") {
			dt$ICD9_string[i] <- paste(dt$ICD9_string[i], paste(ICDs_to_add, collapse=","), sep=",")
		} else {
			dt$ICD9_string[i] <- paste(ICDs_to_add, collapse=",")
		}
	}
	i <- i+1
}

dt$ICD9_string <- gsub("\\,+", "\\,", dt$ICD9_string)
dt$ICD9_string <- gsub("^\\,", "", dt$ICD9_string)
dt$ICD9_string <- gsub("\\,$", "", dt$ICD9_string)

# Filter down to the collection of diseases/disorders that have ICD information possible of being present in UK Biobank

dt_to_lookup <- dt %>% 
	mutate(
		coding_ICD10_string = gsub("\\,", " | ", dt$ICD10_string),
		coding_ICD9_string = gsub("\\,", " | ", dt$ICD9_string)
		) %>% select(Cause, coding_ICD10_string, coding_ICD9_string)
dt_looked_up <- ukb_icd_counts(dt_to_lookup)
dt_looked_up <- data.table(dt_looked_up %>% rename(cause_name = Cause) %>% mutate(UKBB_prevalence = ((`count`/(502413))*100000)))
setkey(dt_looked_up, "cause_name")

dt <- merge(dt, dt_looked_up)
dt <- dt %>% mutate(
	ICD10_string_BRaVa = gsub("([A-Z][0-9]{2})([0-9]+)", "\\1\\.\\2", dt$ICD10_string),
	ICD9_string_BRaVa = gsub("([0-9]{3})([0-9]+)", "\\1\\.\\2", dt$ICD9_string))
dt <- dt %>% mutate(ICD9_string_BRaVa = gsub("(V[0-9]{2})([0-9]+)", "\\1.\\2", dt$ICD9_string_BRaVa))
dt <- dt %>% mutate(ICD_string_BRaVa = paste(ICD10_string_BRaVa, ICD9_string_BRaVa, sep=","))
# Finally, comma separate and combine the ICD9 and ICD10 codes together for inclusion on ICD_phecode tab of BRaVa_Nominate_Phenotypes spreadsheet.
fwrite(dt, file="data/output/GBD_summarised_adult_diseases_with_UKB.tsv", sep='\t')


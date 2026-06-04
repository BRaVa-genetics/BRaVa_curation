# The diseases and their prevalences.
library(xml2)
library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)

a <- read_xml("data/en_product9_prev.xml")
a_list <- as_list(a)

colnames <- names(unlist(a_list$JDBOR$DisorderList[1]))
dt <- setNames(data.table(matrix(unlist(a_list$JDBOR$DisorderList[1]), nrow=1)), colnames)

for (i in 2:length(a_list$JDBOR$DisorderList)) {
	colnames_tmp <- names(unlist(a_list$JDBOR$DisorderList[i]))
	dt <- merge(dt, setNames(data.table(matrix(unlist(a_list$JDBOR$DisorderList[i]), nrow=1)), colnames_tmp))
	print(i)
}

disorder_names <- c()
disorder_names <- c(disorder_names, a_list$JDBOR$DisorderList[i]$Disorder$Name[[1]])
dt_list <- list()
for (i in 1:length(a_list$JDBOR$DisorderList)) {
	dt_list[[i]] <- rbindlist(a_list$JDBOR$DisorderList[i], fill=TRUE)
}

dt_rbind_list <- rbindlist(dt_list)
dt_rbind_list$DisorderGroup <- unlist(dt_rbind_list$DisorderGroup)
dt_rbind_list$DisorderType <- unlist(dt_rbind_list$DisorderType)

dt_rbind_list$point_prevalence <- NA
dt_rbind_list$point_prevalence_class <- NA
dt_rbind_list$location <- NA

for (i in 1:nrow(dt_rbind_list)) {
	print(i)
	if (
		unlist(dt_rbind_list$PrevalenceList[[i]]$PrevalenceType$Name) == "Point prevalence"
		) {
		if (!is.null(unlist(dt_rbind_list$PrevalenceList[[i]]$PrevalenceClass$Name))) {
			dt_rbind_list$point_prevalence_class[i] <- unlist(dt_rbind_list$PrevalenceList[[i]]$PrevalenceClass$Name)
		} else {
			dt_rbind_list$point_prevalence_class[i] <- ""
		}
		dt_rbind_list$point_prevalence[i] <- unlist(dt_rbind_list$PrevalenceList[[i]]$ValMoy)
		dt_rbind_list$location[i] <- unlist(dt_rbind_list$PrevalenceList[[i]]$PrevalenceGeographic$Name)
	}
}

dt_rarish <- dt_rbind_list %>% mutate(point_prevalence = as.numeric(point_prevalence)) %>% filter(point_prevalence >= 4)
dt_rarish <- dt_rarish %>% mutate(Name = unlist(Name), OrphaCode = unlist(OrphaCode), DisorderType = unlist(DisorderType), DisorderGroup = unlist(DisorderGroup))
dt_final_summary <- dt_rarish %>% group_by(Name, OrphaCode, DisorderType, DisorderGroup) %>% summarise(mean_prev = mean(point_prevalence))

# Average age of death for the diseases (probably a subset of them)
b <- read_xml("data/en_product9_ages.xml")
b_list <- as_list(b)

dt <- data.table(
	OrphaCode = unlist(b_list$JDBOR$DisorderList[1]$Disorder$OrphaCode),
	ExpertLink = unlist(b_list$JDBOR$DisorderList[1]$Disorder$ExpertLink),
	Name = unlist(b_list$JDBOR$DisorderList[1]$Disorder$Name),
	DisorderType = unlist(b_list$JDBOR$DisorderList[1]$Disorder$DisorderType$Name),
	DisorderGroup = unlist(b_list$JDBOR$DisorderList[1]$Disorder$DisorderGroup$Name)
	)

if (attr(b_list$JDBOR$DisorderList[1]$Disorder$AverageAgeOfDeathList, "count") != "0") {
	dt$AverageAgeOfDeath <- paste(unlist(b_list$JDBOR$DisorderList[1]$Disorder$AverageAgeOfDeathList), collapse=",")
} else {
	dt$AverageAgeOfDeath <- NA
}

if (attr(b_list$JDBOR$DisorderList[1]$Disorder$AverageAgeOfOnsetList, "count") != "0") {
	dt$AverageAgeOfOnset <- paste(unlist(b_list$JDBOR$DisorderList[1]$Disorder$AverageAgeOfOnsetList), collapse=",")
} else {
	dt$AverageAgeOfOnset <- NA
}

if (attr(b_list$JDBOR$DisorderList[1]$Disorder$TypeOfInheritanceList, "count") != "0") {
	dt$TypeOfInheritance <- paste(unlist(b_list$JDBOR$DisorderList[1]$Disorder$TypeOfInheritanceList), collapse=",")
} else {
	dt$TypeOfInheritance <- NA
}

for (i in 2:length(b_list$JDBOR$DisorderList)) {
	print(i)
	dt_tmp <- data.table(
	OrphaCode = unlist(b_list$JDBOR$DisorderList[i]$Disorder$OrphaCode),
	ExpertLink = unlist(b_list$JDBOR$DisorderList[i]$Disorder$ExpertLink),
	Name = unlist(b_list$JDBOR$DisorderList[i]$Disorder$Name),
	DisorderType = unlist(b_list$JDBOR$DisorderList[i]$Disorder$DisorderType$Name),
	DisorderGroup = unlist(b_list$JDBOR$DisorderList[i]$Disorder$DisorderGroup$Name)
	)
	if (attr(b_list$JDBOR$DisorderList[i]$Disorder$AverageAgeOfDeathList, "count") != "0") {
		dt_tmp$AverageAgeOfDeath <- paste(unlist(b_list$JDBOR$DisorderList[i]$Disorder$AverageAgeOfDeathList), collapse=",")
	} else {
		dt_tmp$AverageAgeOfDeath <- NA
	}
	if (attr(b_list$JDBOR$DisorderList[i]$Disorder$AverageAgeOfOnsetList, "count") != "0") {
		dt_tmp$AverageAgeOfOnset <- paste(unlist(b_list$JDBOR$DisorderList[i]$Disorder$AverageAgeOfOnsetList), collapse=",")
	} else {
		dt_tmp$AverageAgeOfOnset <- NA
	}
	if (attr(b_list$JDBOR$DisorderList[i]$Disorder$TypeOfInheritanceList, "count") != "0") {
		dt_tmp$TypeOfInheritance <- paste(unlist(b_list$JDBOR$DisorderList[i]$Disorder$TypeOfInheritanceList), collapse=",")
	} else {
		dt_tmp$TypeOfInheritance <- NA
	}
	dt <- rbind(dt, dt_tmp)#setNames(data.table(matrix(unlist(b_list$JDBOR$DisorderList[i]), nrow=1)), colnames_tmp))
}

dt_final_summary <- data.table(dt_final_summary)

setkeyv(dt_final_summary, c("Name", "OrphaCode", "DisorderType", "DisorderGroup"))
setkeyv(dt, c("Name", "OrphaCode", "DisorderType", "DisorderGroup"))

dt <- merge(dt, dt_final_summary)

dt <- dt %>% mutate(`Age of death` = dt$AverageAgeOfDeath)
dt$`Age of death`[
	(grepl("child", dt$AverageAgeOfDeath) | 
	grepl("young", dt$AverageAgeOfDeath) |
	grepl("adolescent", dt$AverageAgeOfDeath) |
	grepl("infantile", dt$AverageAgeOfDeath) |
	grepl("stillbirth", dt$AverageAgeOfDeath))
	] <- "Young adult or younger"
dt$`Age of death`[
	is.na(dt$AverageAgeOfDeath) |
	grepl("No data available", dt$AverageAgeOfDeath)
	] <- "No data"
dt$`Age of death`[
	!(grepl("Young adult or younger", dt$AverageAgeOfDeath) |
	  grepl("No data", dt$AverageAgeOfDeath))
	] <- "Adult or normal life expectancy"

# Remove phenotypes with no average age of death information
dt %>% filter(!(DisorderType %in% c(
	"Etiological subtype",
	"Particular clinical situation in a disease or syndrome")
))

pdf(file="plot_all.pdf", width=18, height=25)
p <- ggplot(dt %>% mutate(Name=reorder(reorder(Name, mean_prev), as.numeric(as.factor(DisorderType)))), aes(x=mean_prev, y=Name)) + 
	geom_point(aes(color = DisorderType)) + 
	scale_x_continuous(trans='log10') +
	xlab("Point prevalence (per 100,000)") + 
	ylab("") + 
	theme_classic() + theme(axis.text.y = element_text(size = 6)) +
	geom_vline(xintercept = 25, color="red") + 
	geom_vline(xintercept = 1000, color="red")
print(p)
dev.off()

pdf(file="plot_rareish.pdf", width=10, height=10)
p <- ggplot(dt %>% filter(mean_prev >= 25, mean_prev <= 1000) %>% mutate(Name=reorder(reorder(Name, mean_prev), as.numeric(as.factor(DisorderType)))), aes(x=mean_prev, y=Name)) + 
	geom_point(aes(color = DisorderType)) + 
	scale_x_continuous(trans='log10') +
	xlab("Point prevalence (per 100,000)") + 
	ylab("") + 
	theme_classic() + theme(axis.text.y = element_text(size = 6)) +
	geom_vline(xintercept = 25, color="red") + 
	geom_vline(xintercept = 1000, color="red")
print(p)
dev.off()

fwrite(dt, file="data/output/orphanet_summarised_adult_diseases.tsv", sep='\t')

create_orphanet_list <- function(
	orphanet_icd10_codes="data/icd_mappings/ORPHAnomenclature_MasterFile_en.xlsx",
	skip_orpha=1) {
	dt <- read_excel(orphanet_icd10_codes, skip=skip_orpha) %>% filter(!is.na(ICDcodes))
	# Ensure that the ORPHAcode is now unique
	dt <- dt %>% group_by(ORPHAcode, PreferredTerm) %>% summarise(
		ICDcodes_string = paste(ICDcodes, collapse=" | "))
	dt <- dt %>% mutate(coding_string = gsub("\\.", "", ICDcodes_string))
	return(data.table(dt))
}

# Sanity check that the coding in UKB for ICD10 codes is the same as the actual code, but with the . removed.

ukb_icd10_encoding <- "data/icd_mappings/coding19.tsv"
dt_ukb <- fread(ukb_icd10_encoding) %>% mutate(
	ICDcodes=gsub(" .*", "", meaning),
	coding=gsub("Block ", "", coding)) %>% select(coding, ICDcodes)
dt_ukb[which(gsub("\\.", "", dt_ukb$ICDcodes) != dt_ukb$coding),]
# Yep, all good.

ukb_icd10_counts <- function(
	dt_ICD10_codes_to_lookup,
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
	ICD10s <- c("41202", "41204", "40006", "40001", "40002")
	cols <-  get_cols(ICD10s, dt)
	select_cols <- rep("character", (length(cols) + 1))
	names(select_cols) <- c("eid", cols)

	# Read in the entire file ensuring these columns are encoded as characters to avoid NA weirdness.
	dt <- fread(phenotype_file, na.strings=NULL, select=select_cols)
	ICD10_cols <-  get_cols(ICD10s, dt)

	# Combine vector into a single string
	dt[, ICD10_string := do.call(paste,.SD), .SDcols=ICD10_cols]
	# Remove trailing whitespace
	dt[, ICD10_string := gsub("( )+", " ", ICD10_string)]
	# Pad start and end with space
	dt[, ICD10_string := gsub("^(.*)$", " \\1 ", ICD10_string)]

	dt_ICD10_codes_to_lookup$count <- NA
	for (i in 1:nrow(dt_ICD10_codes_to_lookup)) {
		grep_exp <- dt_ICD10_codes_to_lookup$coding_string[i]
		dt_ICD10_codes_to_lookup$count[i] <- sum(grepl(grep_exp, dt$ICD10_string))
	}
	dt_ICD10_codes_to_lookup <- data.table(dt_ICD10_codes_to_lookup %>% select(OrphaCode, coding_string, count))
	setkey(dt_ICD10_codes_to_lookup, "OrphaCode")
	return(dt_ICD10_codes_to_lookup)
}

dt <- fread("data/output/orphanet_summarised_adult_diseases.tsv", key="OrphaCode")

# Determine the counts of these in the mapped ICD codes, in UK Biobank, and the associated name of the coding.
dt_orpha_icd <- data.table(create_orphanet_list() %>% rename(OrphaCode = ORPHAcode) %>% mutate(OrphaCode = as.integer(OrphaCode)))
setkey(dt_orpha_icd, "OrphaCode")
dt_to_lookup <- merge(dt, dt_orpha_icd)
dt_looked_up <- ukb_icd10_counts(dt_to_lookup)

dt <- merge(dt, dt_looked_up, all.x=TRUE)
dt <- dt %>% mutate(ICD10_string_BRaVa = gsub("([A-Z][0-9]{2})([0-9]+)", "\\1\\.\\2", dt$coding_string))
dt <- dt %>% mutate(ICD10_string_BRaVa = gsub(" \\| ", ",", dt$ICD10_string_BRaVa))

setkey(dt, "Name")
# This file is then pasted into google sheets, and manually curated (sheet 1)
fwrite(dt, file="data/output/orphanet_summarised_adult_diseases_ukb_counts.tsv", sep="\t")



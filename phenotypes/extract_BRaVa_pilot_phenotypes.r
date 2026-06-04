library(googlesheets4)
library(dplyr)
library(data.table)

dt <- read_sheet("https://docs.google.com/spreadsheets/d/1YqdSyxf2OyoIYvLnDVj7NmbpebtppsgyJSq18gkVAWI/edit#gid=1716081249",
	sheet="Sequenced_Sample_Sizes", skip=1)
	
dt <- dt %>% select("...2",
	"...3",
	"...68",
	"...69",
	"...70",
	"...71",
	"...72",
	"UK Biobank...128")
dt$`...68` <- as.character(dt$`...68`)
dt$`...69` <- as.character(dt$`...69`)
dt$`...70` <- as.character(dt$`...70`)
dt$`...71` <- as.character(dt$`...71`)
dt$`...72` <- as.character(dt$`...72`)
names(dt) <- c("Phenotype", "Trait_type", "AFR", "AMR", "EAS", "EUR", "SAS", "Include")

# Clean up names
dt$Phenotype <- gsub(" [ ]+", " ", dt$Phenotype)
dt$Phenotype <- gsub(" / ", "/", dt$Phenotype)
dt$Phenotype <- gsub("([a-z])(\\()", "\\1 \\2", dt$Phenotype)
dt$Phenotype <- gsub(" \\([A-Z,\\/a-z]+\\)", "", dt$Phenotype)
dt$Phenotype <- gsub(",", "", dt$Phenotype)
dt$Phenotype <- gsub("-", " ", dt$Phenotype)
dt$Phenotype <- gsub("[\\(,\\)]", "", dt$Phenotype)
dt$Phenotype <- gsub(" ", "_", dt$Phenotype)

dt <- dt %>% filter(Include == 1)

dt <- dt %>% mutate(
	AFR = as.integer(AFR),
	AMR = as.integer(AMR),
	EAS = as.integer(EAS),
	EUR = as.integer(EUR),
	SAS = as.integer(SAS)
	)

dt_cat <- dt %>% filter(Trait_type != "N/A (quantitative trait)")
dt_cts <- dt %>% filter(Trait_type == "N/A (quantitative trait)")

# Create a string for each ancestry
# Filter case control traits using N case > 100
# Filter continuous traits using N > 500

# Define the set of sex specific phenotypes

remove_from_both_sexes_female <- c(
	"Benign_and_in_situ_cervical_and_uterine_neoplasms",
	"Breast_cancer",
	"Cervical_cancer",
	"Endometriosis",
	"Excess_frequent_and_irregular_menstrual_bleeding",
	"Female_infertility",
	"Female_infertility_anatomic_causes_only",
	"Maternal_hemorrhage",
	"Maternal_hypertensive_disorders",
	"Ovarian_Cancer",
	"Placental_insufficiency",
	"Polycystic_ovarian_syndrome",
	"Preeclampsia",
	"Pregnancy_Loss",
	"Uterine_cancer")

remove_from_both_sexes_male <- c("Male_infertility")

for (anc in c("AFR", "AMR", "EAS", "EUR", "SAS"))
{
	phenos_cat <- unlist((dt_cat %>% filter(dt_cat[[anc]] >= 100))$Phenotype)
	phenos_cts <- unlist((dt_cts %>% filter(dt_cts[[anc]] >= 500))$Phenotype)

	phenos_cat_both <- phenos_cat[!(phenos_cat %in% c(remove_from_both_sexes_female, remove_from_both_sexes_male))]
	phenos_cat_female <- phenos_cat[phenos_cat %in% remove_from_both_sexes_female]
	phenos_cat_male <- phenos_cat[phenos_cat %in% remove_from_both_sexes_male]

	phenos_cts_both <- phenos_cts[!(phenos_cts %in% c(remove_from_both_sexes_female, remove_from_both_sexes_male))]
	phenos_cts_female <- phenos_cts[phenos_cts %in% remove_from_both_sexes_female]
	phenos_cts_male <- phenos_cts[phenos_cts %in% remove_from_both_sexes_male]

	if (length(phenos_cat_both) > 0) {
		fwrite(list(paste(phenos_cat_both, collapse=" ")),
			file=paste0("data/", anc, "_cat_both.txt"), quote=FALSE)
	}
	if (length(phenos_cts_both) > 0) {
		fwrite(list(paste(phenos_cts_both, collapse=" ")),
			file=paste0("data/", anc, "_cts_both.txt"), quote=FALSE)
	}
	if (length(phenos_cat_female) > 0) {
		fwrite(list(paste(phenos_cat_female, collapse=" ")),
			file=paste0("data/", anc, "_cat_female.txt"), quote=FALSE)
	}
	if (length(phenos_cts_female) > 0) {
		fwrite(list(paste(phenos_cts_female, collapse=" ")),
			file=paste0("data/", anc, "_cts_female.txt"), quote=FALSE)
	}
	if (length(phenos_cat_male) > 0) {
		fwrite(list(paste(phenos_cat_male, collapse=" ")),
			file=paste0("data/", anc, "_cat_male.txt"), quote=FALSE)
	}
	if (length(phenos_cts_male) > 0) {
		fwrite(list(paste(phenos_cts_male, collapse=" ")),
			file=paste0("data/", anc, "_cts_male.txt"), quote=FALSE)
	}
}

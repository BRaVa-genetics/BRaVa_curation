library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)

source("meta_analysis_utils.r")

# Short script to compare the number of assocations below threshold detailed in the genebass paper for:

# 1. The ASHG meta-analysis
# 2. The Feb 25 meta-analysis

skat_o_T <- 2.5e-7

# Gene-level results

# ASHG
meta_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/ashg_2024/gene/n_cases_100/"
meta_files <- dir(meta_dir, pattern = ".gz$")
meta_files <- meta_files[-grep("cauchy", meta_files)]

dt_list <- list()
for (file in paste0(meta_dir, "/", meta_files)) {
	dt_list[[file]] <- fread(file) %>% 
		filter(max_MAF != "0.01", Pvalue < 2.5e-7) %>%
		filter(Group %in% c(
			"damaging_missense_or_protein_altering",
			"pLoF",
			"pLoF;damaging_missense_or_protein_altering"
			)) %>% mutate(filename=file)
}
dt_ashg <- rbindlist(dt_list)
dt_ashg$freeze <- "ashg_2024"

# Feb 25
meta_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/"
meta_files <- dir(meta_dir, pattern = ".gz$")

dt_list <- list()
for (file in paste0(meta_dir, "/", meta_files)) {
	dt_list[[file]] <- fread(file) %>% 
		filter(max_MAF != "0.01", Pvalue < 2.5e-7) %>%
		filter(Group %in% c(
			"damaging_missense_or_protein_altering",
			"pLoF",
			"pLoF;damaging_missense_or_protein_altering"
			)) %>% mutate(filename=file)
}
dt_feb25 <- rbindlist(dt_list)
dt_feb25$freeze <- "feb_2025"

dt <- rbind(dt_ashg, dt_feb25)
dt <- dt %>% mutate(phenotype = gsub(".*\\/([A-Z0-9a-z]*)_.*", "\\1", filename))

fwrite(dt, file="/well/lindgren/dpalmer/BRaVa_curation/meta_analysis/summary_of_significant_genes.tsv.gz", sep="\t", quote=FALSE)
dt <- fread("/well/lindgren/dpalmer/BRaVa_curation/meta_analysis/summary_of_significant_genes.tsv.gz")
# Our rule is going to be Inverse variance weighted if we can (burden), otherwise Stouffer

# Stouffer
Stouffer_hits <- unique(dt %>%
	filter(type == "Stouffer", class != "Burden") %>%
	select(Region, phenotype, freeze) %>%
	filter(!phenotype %in% c("MatHem", "AlcCons", "ALT", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI"),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))))
Stouffer_hits %>% group_by(freeze) %>% summarise(count = n())
print(Stouffer_hits %>% group_by(phenotype, freeze) %>% summarise(count = n()), n=100)

Fisher_hits <- unique(dt %>%
	filter(type == "Weighted Fisher") %>%
	select(Region, phenotype, freeze) %>%
	filter(!phenotype %in% c("MatHem", "AlcCons", "ALT", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI"),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))))
Fisher_hits %>% group_by(freeze) %>% summarise(count = n())
print(Fisher_hits %>% group_by(phenotype, freeze) %>% summarise(count = n()), n=100)

IVW_hits <- unique(dt %>%
	filter(type == "Inverse variance weighted") %>%
	select(Region, phenotype, freeze) %>%
	filter(!phenotype %in% c("MatHem", "AlcCons", "ALT", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI"),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))))
IVW_hits %>% group_by(freeze) %>% summarise(count = n())
print(IVW_hits %>% group_by(phenotype, freeze) %>% summarise(count = n()), n=100)

Stouffer_and_IVW_hits <- unique(rbind(IVW_hits, Stouffer_hits))
Stouffer_and_IVW_hits %>% group_by(freeze) %>% summarise(count = n())

setdiff((Stouffer_and_IVW_hits %>%
	filter(freeze == "feb_2025") %>%
	mutate(hit = paste(phenotype, Region)))$hit,
	(Stouffer_and_IVW_hits %>%
	filter(freeze == "ashg_2024") %>%
	mutate(hit = paste(phenotype, Region)))$hit
	)

All_hits <- unique(dt %>%
	select(Region, phenotype, freeze) %>%
	filter(!phenotype %in% c("MatHem", "AlcCons", "ALT", "AST", "BMI", "HDLC", "Height", "LDLC", "TChol", "TG", "WHRBMI"),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))))
All_hits %>% group_by(freeze) %>% summarise(count = n())
print(All_hits %>% group_by(phenotype, freeze) %>% summarise(count = n()), n=100)







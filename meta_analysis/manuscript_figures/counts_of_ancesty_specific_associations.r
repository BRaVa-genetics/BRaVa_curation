# The following is the ancestry specific equivalent of the above.

# 197 continuous and 90 binary - this is the local count.
# Next, login and run this code on each of the genetic ancestry labels.
# Location: /well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100
files_list <- c("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AMR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/AFR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/EUR",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/SAS",
	"/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100/non_EUR")

meta_list <- list()
for (file in files_list)
{
	files <- dir(file, full.names=TRUE)
	files <- grep("tsv.gz", files, value=TRUE)
	files <- files[-which(grepl("extra_cauchy", files))]
	print(files)
	meta_list[[file]] <- list()

	for (f in files) {
		cat(f, "\n")
		phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", f)
		meta_list[[file]][[f]] <- fread(f) %>% filter(
			max_MAF %in% c("1e-04", "0.001"),
			Group %in% c(
				"damaging_missense_or_protein_altering",
				"pLoF",
				"pLoF;damaging_missense_or_protein_altering")) %>%
		mutate(phenotype = phenotype) %>% filter(
			(Pvalue < 6.7e-7 & class == "Burden" & type == "Inverse variance weighted") |
			(Pvalue < 2.5e-7 & class %in% c("SKAT", "SKAT-O") & type == "Stouffer"))
	}

	meta_list[[file]] <- rbindlist(meta_list[[file]]) %>% mutate(case_control =
		ifelse(phenotype %in% case_ctrl, TRUE,
			ifelse(phenotype %in% cts, FALSE, NA)))
	if (nrow(meta_list[[file]]) == 0) {
		cat("No significant associations\n")
		next
	}
	meta_list_unique <- meta_list[[file]] %>% group_by(case_control) %>% 
		filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
		summarise(count = length(unique(paste(Region, phenotype))))
	print(meta_list_unique)
}

for (anc in names(meta_list)) {
	meta_list[[anc]]$ancestry <- ifelse(anc == "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
		"all", gsub(".*/", "", anc))
}
meta <- rbindlist(meta_list)
# Now...
for (anc in c("AMR", "AFR", "EAS", "EUR", "non_EUR", "SAS", "all")) {
	print(
		setdiff(
			meta %>% filter(ancestry == !!anc) %>% select(phenotype, Region),
			meta %>% filter(ancestry == "EUR") %>% select(phenotype, Region))
	)
}

# This is contrasting where variance based testing hits things that are missed by burden testing.
setdiff(meta %>% filter(ancestry == "all") %>% 
	filter(Group %in% c("pLoF;damaging_missense_or_protein_altering",
		"damaging_missense_or_protein_altering")) %>% select(phenotype, Region),
	meta %>% filter(ancestry == "all") %>% 
	filter(Group == "pLoF") %>% select(phenotype, Region))
# 68, 77 height

setdiff(meta %>% filter(ancestry == "all") %>% 
	filter(Group %in% c("damaging_missense_or_protein_altering")) %>% select(phenotype, Region),
	meta %>% filter(ancestry == "all") %>% 
	filter(Group == "pLoF") %>% select(phenotype, Region))
# 28, 32 height

setdiff(meta %>% filter(ancestry == "all") %>% 
	filter(class %in% c("SKAT", "SKAT-O")) %>% select(phenotype, Region),
	meta %>% filter(ancestry == "all") %>% 
	filter(class == "Burden") %>% select(phenotype, Region))
# 19, 20 height
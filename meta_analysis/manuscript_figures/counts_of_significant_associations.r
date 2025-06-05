library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)

source("meta_analysis_utils.r")

# Short script to compare the number of assocations below threshold detailed in the genebass paper for:

# 1. UK-Biobank alone
# 2. The meta-analysis up to yesterday
# 3. The meta-analysis up to today

skat_o_T <- 2.5e-7

# Gene-level results

# UK-Biobank
ukb_dir <- "../data/meta_analysis/gcloud/uk-biobank/cleaned/gene"
files <- dir(ukb_dir)
dt_list <- list()
for (f in files) {
	file_info <- extract_file_info(f)
	if (!file_info$binary | file_info$ancestry != "EUR") {
		next
	}
	cat(paste0(file_info$phenotype, "...", file_info$ancestry, "\n"))
	dt_tmp <- fread(paste0(ukb_dir, "/", f))
	dt_tmp <- dt_tmp %>% mutate(phenotypeID = file_info$phenotype)
	dt_list[[f]] <- dt_tmp %>% 
		filter(
			Group %in% c("pLoF", "Cauchy"),
			((max_MAF == 0.001) | (is.na(max_MAF)))
		) %>% 
		filter(
			(Pvalue < skat_o_T) |
			(Pvalue_Burden < skat_o_T) | 
			(Pvalue_SKAT < skat_o_T)
		) %>% 
		select(Region, Group, max_MAF,
			Pvalue, Pvalue_Burden, Pvalue_SKAT,
			phenotypeID)
}

dt_ukb <- rbindlist(dt_list)

dt_ukb_count <- dt_ukb %>% group_by(Group) %>%
	summarise(skat_o_count = sum(Pvalue < skat_o_T),
		skat_count = sum(Pvalue_SKAT < skat_o_T),
		burden_count = sum(Pvalue_Burden < skat_o_T))

# Meta-analysis results (v3)
meta_dir <- "../data/meta_analysis/meta_results/n_cases_100"
files <- grep("gz$", dir(meta_dir), value=TRUE)
files <- files[-which(grepl("extra", files))]
dt_list <- list()
for (f in files) {
	phenotype <- gsub("_.*", "", f)
	if (phenotype %in% c(
		"ALT", "AlcCons", "AST", "BMI", "CRP", "HDLC", "Height",
		"LDLC", "TChol", "TG", "WHRBMI", "MatHem")) {
		next
	}
	cat(paste0(phenotype, "\n"))
	dt_tmp <- fread(paste0(meta_dir, "/", f))
	dt_tmp <- dt_tmp %>% mutate(phenotypeID = phenotype)
	dt_list[[f]] <- dt_tmp %>% 
		filter(
			Group %in% c("pLoF", "Cauchy"),
			((max_MAF == "0.001") | (max_MAF == "Cauchy")),
			type == "Weighted Fisher",
			Pvalue < skat_o_T
		) %>% select(
		Region, Group, max_MAF, Pvalue, class,
		phenotypeID
		)
}

dt_meta_v3 <- rbindlist(dt_list)
dt_meta_v3_count <- dt_meta_v3 %>% group_by(Group, class) %>% 
	summarise(count = sum(Pvalue < skat_o_T))


# Plotting power gain starts here
ukb_dir <- "../data/meta_analysis/gcloud/uk-biobank/cleaned/gene"
ukb_files <- dir(ukb_dir)[grep("\\.EUR\\.", dir(ukb_dir))]
meta_dir <- "../data/meta_analysis/meta_results/n_cases_100"
meta_files <- dir(meta_dir, pattern = ".gz$")
meta_files <- meta_files[-grep("cauchy", meta_files)]

# Determine the overall effective N weighting by sqrt of effective sample size of each of the studies
files <- grep("cleaned.*gene.*", dir("../data/meta_analysis/gcloud/", recursive=TRUE), value=TRUE)
information <- sapply(files, extract_file_info)
information <- rbindlist(information, fill=TRUE)
information_ukb <- information %>% filter(
	dataset == "uk-biobank/gene/uk-biobank",
	ancestry=="EUR"
	)
get_n_eff <- function(dt) {
	return(
		rbind(
			dt %>% filter(binary) %>% group_by(phenotype, binary) %>% 
			mutate(n_eff = 4/((1/as.integer(n_cases)) + (1/as.integer(n_controls)))) %>%
			summarise(meta_n_eff = sum(n_eff)),
			dt %>% filter(!binary) %>% group_by(phenotype, binary) %>%
			mutate(n_eff = as.integer(n)) %>% summarise(meta_n_eff = sum(n_eff))
			)
	)
}

dt_neff_ukb <- data.table(get_n_eff(information_ukb), key=c("phenotype", "binary"))[, n_eff_ukb := meta_n_eff]
dt_neff_ukb[, meta_n_eff := NULL]
dt_neff_meta <- data.table(get_n_eff(information), key=c("phenotype", "binary"))[, n_eff_meta := meta_n_eff]
dt_neff_meta[, meta_n_eff := NULL]
dt_neff <- merge(dt_neff_ukb, dt_neff_meta)
dt_neff[, scaling := n_eff_meta/n_eff_ukb]

dt_list <- list()
for (i in 1:length(ukb_files))
{
	info <- extract_file_info(ukb_files[i])
	meta_file <- grep(paste0("^", info$phenotype, "_"), meta_files, value=TRUE)
	cat(paste(meta_file, collapse="\n"))
	if (length(meta_file) == 0) { 
		print(paste0("skipping ", info$phenotype))
		next
	}
	# Contrast the power between the meta-analysis and the UK Biobank version, and split by phenotype and the use the n_eff
	dt_meta <- fread(paste0(meta_dir, "/", meta_file)) %>%
		filter(max_MAF %in% c("1e-04", "0.001"))
	dt_meta[, max_MAF := as.numeric(max_MAF)]
	dt_ukb <- fread(paste0(ukb_dir, "/", ukb_files[i])) %>%
		filter(max_MAF %in% c("1e-04", "0.001"))
	dt_ukb <- melt(
	  dt_ukb,
	  id.vars = c("Region", "Group", "max_MAF"),
	  measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
	  variable.name = "class",
	  value.name = "Pvalue"
	)

	dt_ukb[class == "Pvalue", class := "SKAT-O"]
	dt_ukb[class == "Pvalue_Burden", class := "Burden"]
	dt_ukb[class == "Pvalue_SKAT", class := "SKAT"]
	dt_ukb <- dt_ukb %>% 
		filter(Group == "pLoF")
	dt_ukb[, Pvalue_ukb := Pvalue]
	dt_ukb[, max_MAF := as.numeric(max_MAF)]
	dt_ukb[, chi2_ukb := qchisq(Pvalue, 1, lower=FALSE)]
	dt_ukb[, Pvalue := NULL]

	setkeyv(dt_ukb, c("Region", "Group", "max_MAF", "class"))
	setkeyv(dt_meta, c("Region", "Group", "max_MAF", "class", "type"))

	dt_meta[, chi2_meta := qchisq(Pvalue, 1, lower=FALSE)]
	dt_list[[info$phenotype]] <- merge(dt_ukb, dt_meta)
	dt_list[[info$phenotype]]$phenotype <- info$phenotype
	dt_list[[info$phenotype]]$n_eff <- ifelse(info$binary, 4/(1/as.integer(info$n_cases) + 1/as.integer(info$n_controls)), as.integer(info$n))
}

# Stouffer, burden
dt_plot <- rbindlist(dt_list) %>% 
	# filter(Pvalue < 2.5e-7) %>%
	filter(`Pvalue_ukb` < 2.5e-7 | `Pvalue` < 2.5e-7) %>%
	filter(Group == "pLoF", class=="Burden", type=="Inverse variance weighted")
setkey(dt_neff, "phenotype")

renaming_phenotype_list[["AlcCons"]] <- "Alcohol consumption\n(drinks per week)"
renaming_phenotype_list[["AST"]] <- "Aspartate\naminotransferase"
renaming_phenotype_list[["BenIntNeo"]] <- "Benign and in situ\nintestinal neoplasms"
renaming_phenotype_list[["COPD"]] <- "Chronic obstructive\npulmonary disease"
renaming_phenotype_list[["ILDSarc"]] <- "Interstitial lung disease\nand pulmonary sarcoidosis"
renaming_phenotype_list[["EFRMB"]] <- "Excess, frequent and irregular\nmenstrual bleeding"
renaming_phenotype_list[["NonRheuValv"]] <- "Non-rheumatic valvular\nheart disease"
renaming_phenotype_list[["WHRBMI"]] <- "Waist to hip ratio\nadjusted for BMI"

dt_plot <- dt_plot %>% mutate(phenotype = unlist(renaming_phenotype_list[phenotype]))
dt_neff <- dt_neff %>% mutate(phenotype = unlist(renaming_phenotype_list[phenotype]))

# All of them
p <- ggplot(dt_plot, aes(x=chi2_ukb, y=chi2_meta)) + geom_point() + expand_limits(x = 0, y = 0) + facet_wrap(~phenotype, scales='free')
p <- p + geom_abline(data=dt_neff %>% filter(phenotype %in% unique(dt_plot$phenotype)), aes(intercept=0, slope=scaling), col='red') +
geom_abline(data=dt_neff %>% filter(phenotype %in% unique(dt_plot$phenotype)), aes(intercept=0, slope=1), col='grey40', linetype='dashed')
p <- p + theme_minimal() + labs(
	x = TeX("$\\chi^2_{UK-Biobank,(EUR)}"),
	y = TeX("$\\chi^2_{Meta}"))
# pdf(file="chi2_examples.pdf", width=9, height=3)
print(p)
dev.off()

# Nicely behaved examples
dt_plot_examples <- dt_plot %>% filter(phenotype %in% c("Atrial fibrillation", "Breast cancer", "Total cholesterol"))

p <- ggplot(dt_plot_examples, aes(x=chi2_ukb, y=chi2_meta)) + geom_point() + expand_limits(x = 0, y = 0) + facet_wrap(~phenotype, scales='free')
p <- p + geom_abline(data=dt_neff %>% filter(phenotype %in% unique(dt_plot_examples$phenotype)), aes(intercept=0, slope=scaling), col='red') +
geom_abline(data=dt_neff %>% filter(phenotype %in% unique(dt_plot_examples$phenotype)), aes(intercept=0, slope=1), col='grey40', linetype='dashed')
p <- p + theme_minimal() + labs(
	x = TeX("$\\chi^2_{UK-Biobank,(EUR)}"),
	y = TeX("$\\chi^2_{Meta}"))
pdf(file="chi2_examples.pdf", width=7, height=2.5)
print(p)
dev.off()

	# geom_abline(data=dt_neff %>% filter(phenotype %in% unique(dt_plot$phenotype)), aes(intercept=0, slope=1, col='red'))


# For ASHG abstract
dt <- rbindlist(dt_list) %>%
	filter((`Pvalue_ukb` < 1e-10), ((`Pvalue_ukb` > 0) & (`Pvalue` > 0))) %>%
	filter(Group == "pLoF", type=="Stouffer", class=="Burden")
dt <- dt %>% filter(phenotype != "AlcCons", !(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

dt_summary <- merge(dt, dt_neff) %>% mutate(
	scaling_chi = chi2_meta/chi2_ukb,
	implied_increase = ((scaling_chi * n_eff_ukb) - n_eff_ukb),
	actual_increase = n_eff_meta - n_eff_ukb
	) %>% 
	group_by(phenotype, binary) %>% summarise(
		change = median(scaling_chi/scaling),
		scaling_chi_median = median(scaling_chi),
		scaling_median = median(scaling),
		implied_increase_median = median(implied_increase),
		actual_increase_median = median(actual_increase)
		)


dt <- rbindlist(dt_list) %>%
	filter((`Pvalue_ukb` < 1e-10), ((`Pvalue_ukb` > 0) & (`Pvalue` > 0))) %>%
	filter(Group == "pLoF")
dt <- dt %>% filter(phenotype != "AlcCons", !(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

dt_summary <- merge(dt, dt_neff) %>% mutate(scaling_chi = chi2_meta/chi2_ukb) %>% 
	group_by(phenotype, binary) %>% summarise(
		change = median(scaling_chi/scaling),
		scaling_chi = median(scaling_chi),
		scaling = median(scaling)
		)

median((dt_summary %>% filter(binary))$scaling_chi)
median((dt_summary %>% filter(binary))$scaling)

mean((dt_summary %>% filter(binary))$scaling_chi)
mean((dt_summary %>% filter(binary))$scaling)


meta_dir <- "../data/meta_analysis/meta_results/n_cases_100"
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

# Stouffer
unique(rbindlist(dt_list) %>% filter(type == "Stouffer") %>% select(Region, filename)) %>% filter(
	!grepl("AlcCons", filename),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

# Fisher
unique(rbindlist(dt_list) %>% filter(type == "Weighted Fisher") %>% select(Region, filename)) %>% filter(
	!grepl("AlcCons", filename),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

# Inverse variance weighted
unique(rbindlist(dt_list) %>% filter(type == "Inverse variance weighted") %>% select(Region, filename)) %>% filter(
	!grepl("AlcCons", filename),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

# All of them 
unique(rbindlist(dt_list) %>% select(Region, filename)) %>% filter(
	!grepl("AlcCons", filename),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))

dt_check <- rbindlist(dt_list)
dt_het <- dt_check %>% filter(
	!grepl("AlcCons", filename),
	!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))
	) %>% filter(Pvalue_het < 2.5e-7)





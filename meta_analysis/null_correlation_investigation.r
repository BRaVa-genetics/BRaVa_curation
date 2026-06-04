# Investigation into the CAC p-values for synonymous variants

# Loop over phenotypes - extract the phenotypes that were analysed by looking into the
# meta-analysis folder
# Grab all the biobanks
# Extract the p-values, split by MAF, synonymous
source("meta_analysis_utils.r")
library(dplyr)
library(data.table)

sumstats_files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", recursive=TRUE, full.names=TRUE)
sumstats_files <- grep("cleaned.*\\.gene\\.", sumstats_files, value=TRUE)
dt_sumstats <- rbindlist(lapply(sumstats_files, extract_file_info),fill=TRUE)
dt_sumstats <- dt_sumstats %>% dplyr::mutate(dataset = gsub(".*/", "", dataset))
dt_sumstats$file <- sumstats_files
phenotypes <- unique(dt_sumstats$phenotype)

cor_and_count_list <- list()

# Go through all of the largest cohorts for each genetic ancestry, and determine the CAC
# Compare to true CAC for those that we have (validates the approach)
dt_sumstats_binary <- dt_sumstats %>% filter(binary, software != "regenie")
dt_sumstats_binary_summary <- dt_sumstats_binary %>%
  group_by(ancestry) %>%
  summarise(
    n_max = max(as.integer(n_cases)+as.integer(n_controls)),
    dataset_n = dataset[which.max(as.integer(n_cases)+as.integer(n_controls))],
    file_n = file[which.max(as.integer(n_cases)+as.integer(n_controls))],
    .groups = "drop"
  )

dt_sumstats_cts <- dt_sumstats %>% filter(!binary, software != "regenie")
dt_sumstats_cts_summary <- dt_sumstats_cts %>%
  group_by(ancestry) %>%
  summarise(
    n_max = max(as.integer(n), na.rm = TRUE),
    dataset_n = dataset[which.max(as.integer(n))],
    file_n = file[which.max(as.integer(n))],
    .groups = "drop"
  )

# Determine the collection of samples with the largest sample size for each ancestry to use as 
# a reference for those where MAC is not available.
# Let's determine references

# Let's contrast to AoU (MAC >= 20) as a comparison

# AFR reference - AoU (created already)
# AMR reference - AoU (created already)
# Note - this was done as I incorrectly set the CAF - half it. Don't rerun this!
# aou_files <- grep("/reference.*.all-of-us*", dir("manuscript_tables", full.names=TRUE), value=TRUE)
# for (f in aou_files) {
# 	dt <- fread(f) %>% mutate(CAF = CAF/2)
# 	fwrite(dt, file=f, sep='\t')
# }

# EAS reference - BBJ or AoU
info <- dt_sumstats_binary_summary %>% filter(ancestry == "EAS")
n_max <- info$n_max
file <- info$file_n
dataset <- info$dataset_n
dt_binary_ref <- fread(file) %>%
	filter(max_MAF == 0.001, Group=="synonymous") %>% 
		mutate(CAF=MAC/(2*n_max)) %>% 
		select(Region, Group, max_MAF, MAC, CAF) %>% rename(MAC_ref=MAC)
setkeyv(dt_binary_ref, c("Region", "Group", "max_MAF"))
fwrite(dt_binary_ref, file = paste0("manuscript_tables/reference.EAS.", n_max, ".", dataset, ".tsv.gz"), sep='\t')

# EUR reference - UKB
info <- dt_sumstats_cts_summary %>% filter(ancestry == "EUR")
n_max <- info$n_max
file <- info$file_n
dataset <- info$dataset_n
dt_cts_ref <- fread(file) %>% 
	filter(max_MAF == 0.001, Group=="synonymous") %>% 
	mutate(CAF=MAC/(2*n_max)) %>% 
	select(Region, Group, max_MAF, MAC, CAF) %>% rename(MAC_ref=MAC)
setkeyv(dt_cts_ref, c("Region", "Group", "max_MAF"))
fwrite(dt_cts_ref, file = paste0("manuscript_tables/reference.EUR.", n_max, ".", dataset, ".tsv.gz"), sep='\t')

# SAS reference G&H
info <- dt_sumstats_cts_summary %>% filter(ancestry == "SAS")
n_max <- info$n_max
file <- info$file_n
dataset <- info$dataset_n
dt_cts_ref <- fread(file) %>% 
	filter(max_MAF == 0.001, Group=="synonymous") %>% 
	mutate(CAF=MAC/(2*n_max)) %>% 
	select(Region, Group, max_MAF, MAC, CAF) %>% rename(MAC_ref=MAC)
setkeyv(dt_cts_ref, c("Region", "Group", "max_MAF"))
fwrite(dt_cts_ref, file = paste0("manuscript_tables/reference.SAS.", n_max, ".", dataset, ".tsv.gz"), sep='\t')

# Let's contrast for UKB (MAC >= 20) as a comparison
dt_ref_ukb_eur <- fread("manuscript_tables/reference.EUR.402167.uk-biobank.tsv.gz")
dt_ref_aou_eur <- fread("manuscript_tables/reference.EUR.210399.all-of-us.tsv.gz")

setkeyv(dt_ref_ukb_eur, c("Region", "Group", "max_MAF"))
setkeyv(dt_ref_aou_eur, c("Region", "Group", "max_MAF"))
dt_eur <- merge(dt_ref_aou_eur, dt_ref_ukb_eur)
plot(dt_eur$CAF.x, dt_eur$CAF.y)
plot(log10(dt_eur$CAF.x), log10(dt_eur$CAF.y))
abline(0,1, col='red')

# Now EAS BBJ vs EAS AoU
dt_ref_bbj_eas <- fread("manuscript_tables/reference.EAS.10245.bbj.tsv.gz")                     
dt_ref_aou_eas <- fread("manuscript_tables/reference.EAS.8895.all-of-us.tsv.gz")

setkeyv(dt_ref_bbj_eas, c("Region", "Group", "max_MAF"))
setkeyv(dt_ref_aou_eas, c("Region", "Group", "max_MAF"))
dt_eas <- merge(dt_ref_aou_eas, dt_ref_bbj_eas)
plot(dt_eas$CAF.x, dt_eas$CAF.y)
plot(log10(dt_eas$CAF.x), log10(dt_eas$CAF.y))
abline(0,1, col='red')

# Now SAS G&H vs SAS AoU
dt_ref_gh_sas <- fread("manuscript_tables/reference.SAS.32043.genes-and-health.tsv.gz")                     
dt_ref_aou_sas <- fread("manuscript_tables/reference.SAS.3794.all-of-us.tsv.gz")

setkeyv(dt_ref_gh_sas, c("Region", "Group", "max_MAF"))
setkeyv(dt_ref_aou_sas, c("Region", "Group", "max_MAF"))
dt_sas <- merge(dt_ref_aou_sas, dt_ref_gh_sas)
plot(dt_sas$CAF.x, dt_sas$CAF.y)
plot(log10(dt_sas$CAF.x), log10(dt_sas$CAF.y))
abline(0,1, col='red')

# Create an overall reference file
# label ancestry, sample size, dataset - then use this to do prediction
ref_files <- grep("reference\\.", dir("manuscript_tables", full.names=TRUE), value=TRUE)
dt_list <- list()
for (file in ref_files) {
	dt_list[[file]] <- fread(file)
	dt_list[[file]]$ancestry <- gsub(".*\\.([A-Z]+)\\..*", "\\1", file)
	dt_list[[file]]$dataset <- gsub(".*\\.([a-z-]+)\\.tsv\\.gz", "\\1", file)
	dt_list[[file]]$n_ref <- as.integer(gsub(".*\\.([0-9]+)\\..*", "\\1", file))
}
dt_ref <- rbindlist(dt_list, use.names=TRUE)
fwrite(dt_ref, "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/reference_cafs.tsv.gz", sep='\t')
dt_ref <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/reference_cafs.tsv.gz")
dt_ref <- dt_ref[, .SD[n_ref == max(as.integer(n_ref))], by = ancestry]
setkeyv(dt_ref, c("Region", "Group", "max_MAF", "ancestry"))

for (i in 1:nrow(dt_sumstats_cts)) {
	dt_tmp <- fread(dt_sumstats_cts$file[i])
	if (!("MAC" %in% names(dt_tmp)))
		next
	anc_tmp <- dt_sumstats_cts$ancestry[i]
	dt_tmp <- dt_tmp %>% filter(Group == "synonymous", max_MAF == 0.001) %>%
		select(Region, Group, max_MAF, MAC) %>%
		mutate(n=as.integer(dt_sumstats_cts$n[i]), ancestry=anc_tmp)
	setkeyv(dt_tmp, c("Region", "Group", "max_MAF", "ancestry"))
	dt_tmp <- merge(dt_tmp, dt_ref)
	if (nrow(dt_tmp) < 100) { next }
	dt_tmp <- dt_tmp %>% mutate(E_MAC = 2*n*CAF)
	print(paste0(dt_sumstats_cts$dataset[i], ", ", dt_sumstats_cts$ancestry[i], ", ", dt_sumstats_cts$phenotype[i]))
	print(cor(dt_tmp$E_MAC, dt_tmp$MAC))
	print(summary(dt_tmp$E_MAC/dt_tmp$MAC))
}

# Now do the same for binary traits
for (i in 1:nrow(dt_sumstats_binary)) {
	dt_tmp <- fread(dt_sumstats_binary$file[i])
	if (!("MAC_case" %in% names(dt_tmp)))
		next
	anc_tmp <- dt_sumstats_binary$ancestry[i]
	dt_tmp <- dt_tmp %>% filter(Group == "synonymous", max_MAF == 0.001) %>%
		select(Region, Group, max_MAF, MAC, MAC_case, MAC_control) %>%
		mutate(n_cases=as.integer(dt_sumstats_binary$n_cases[i]),
			n_controls=as.integer(dt_sumstats_binary$n_controls[i]),
			ancestry=anc_tmp)
	setkeyv(dt_tmp, c("Region", "Group", "max_MAF", "ancestry"))
	dt_tmp <- merge(dt_tmp, dt_ref)
	if (nrow(dt_tmp) < 100) { next }
	dt_tmp <- dt_tmp %>% mutate(
		E_MAC = 2*(n_cases+n_controls)*CAF,
		E_MAC_control = 2*n_controls*CAF,
		E_MAC_case = 2*n_cases*CAF)
	print(paste0(dt_sumstats_binary$dataset[i], ", ", dt_sumstats_binary$ancestry[i], ", ", dt_sumstats_binary$phenotype[i]))
	print(cor(dt_tmp$E_MAC, dt_tmp$MAC))
	print(cor(dt_tmp$E_MAC_control, dt_tmp$MAC_control))
	print(cor(dt_tmp$E_MAC_case, dt_tmp$MAC_case))
	print(summary(dt_tmp$E_MAC/dt_tmp$MAC))
}

# Then use that to paste in the estimated CAC which we can then filter on
dt_cor_count <- data.table(file_1 = character(), file_2 = character(),
	cor = numeric(), pvalue = numeric(), count=integer())
for (pheno in phenotypes) {
	cat(paste0(pheno, "..."))
	for (anc in c("AFR", "AMR", "EAS", "EUR", "SAS")) {
		cat(paste0(anc, "..."))
		dt_tmp <- dt_sumstats %>% filter(ancestry == anc, phenotype == pheno)
		if (nrow(dt_tmp) == 0) {next}
		# Now, loop over the files and determine if the CAC is available
		dt_list <- list()
		if (dt_tmp$binary[1]) {
			for (i in 1:nrow(dt_tmp)) {
				dt_list[[dt_tmp$file[i]]] <- fread(dt_tmp$file[i]) %>% 
					filter(Group == "synonymous", max_MAF == 1e-3)
				dt_list[[dt_tmp$file[i]]]$dataset <- dt_tmp$dataset[i]
				dt_list[[dt_tmp$file[i]]]$phenotype <- pheno
				if ("MAC_case" %in% names(dt_list[[dt_tmp$file[i]]])) {
					dt_list[[dt_tmp$file[i]]] <- dt_list[[dt_tmp$file[i]]] %>% 
						filter(MAC_case > 30, MAC_case < 1000)
				} else {
					dt_list[[dt_tmp$file[i]]] <- dt_list[[dt_tmp$file[i]]] %>%
						mutate(n_cases=as.integer(dt_tmp$n_cases[i]),
							n_controls=as.integer(dt_tmp$n_controls[i]),
							ancestry = dt_tmp$ancestry[i])
					setkeyv(dt_list[[dt_tmp$file[i]]], c("Region", "Group", "max_MAF", "ancestry"))
					dt_list[[dt_tmp$file[i]]] <- merge(dt_list[[dt_tmp$file[i]]], dt_ref) %>% 
						mutate(MAC = 2*(n_cases+n_controls)*CAF,
							MAC_case= 2*n_cases*CAF,
							MAC_control = 2*n_controls*CAF) %>%
						filter(MAC_case > 30, MAC_case < 1000)
				}
				setkeyv(dt_list[[dt_tmp$file[i]]], c("max_MAF", "Region"))
			}
		} else {
			for (i in 1:nrow(dt_tmp)) {
				dt_list[[dt_tmp$file[i]]] <- fread(dt_tmp$file[i]) %>% 
					filter(Group == "synonymous", max_MAF == 1e-3)
				dt_list[[dt_tmp$file[i]]]$dataset <- dt_tmp$dataset[i]
				dt_list[[dt_tmp$file[i]]]$phenotype <- pheno
				if ("MAC" %in% names(dt_list[[dt_tmp$file[i]]])) {
					dt_list[[dt_tmp$file[i]]] <- dt_list[[dt_tmp$file[i]]] %>% 
						filter(MAC > 50)
				} else {
					dt_list[[dt_tmp$file[i]]] <- dt_list[[dt_tmp$file[i]]] %>%
						mutate(n=as.integer(dt_tmp$n[i]),
							ancestry = dt_tmp$ancestry[i])
					setkeyv(dt_list[[dt_tmp$file[i]]], c("Region", "Group", "max_MAF", "ancestry"))
					dt_list[[dt_tmp$file[i]]] <- merge(dt_list[[dt_tmp$file[i]]], dt_ref) %>% 
						mutate(MAC = 2*n*CAF) %>% filter(MAC > 50)
				}
				setkeyv(dt_list[[dt_tmp$file[i]]], c("max_MAF", "Region"))
			}
		}

		if (length(dt_list) == 1) { next }

		for (i in 1:(length(dt_list)-1)) {
			for (j in (i+1):length(dt_list)) {
				dt_ij <- merge(dt_list[[i]], dt_list[[j]])
				if (nrow(dt_ij) < 100) { next } # Fewer than 100 genes to compare
				Z_i <- qnorm(dt_ij$`Pvalue.x`/2) * sign(dt_ij$`BETA_Burden.x`)
				Z_j <- qnorm(dt_ij$`Pvalue.y`/2) * sign(dt_ij$`BETA_Burden.y`)
				test <- cor.test(Z_i, Z_j, alternative = "greater")
				dt_tmp <- data.table(file_1 = names(dt_list)[i], file_2 = names(dt_list)[j],
						cor = test$estimate, pvalue = test$`p.value`, count = nrow(dt_ij))
				dt_cor_count <- rbind(dt_cor_count, dt_tmp)
				print(dt_tmp)
			}
		}
	}
}

# I can use this to generate correlation matrices to pass, or evaluate them on the fly.
dt_sumstats <- dt_sumstats %>% mutate(old=ifelse(grepl("cleaned_v1", file), TRUE, FALSE))
dt_cor_count <- dt_cor_count %>%
  left_join(dt_sumstats %>% rename_with(~ paste0(.x, "_1"), -file),
    by = c("file_1" = "file")) %>%
  left_join(dt_sumstats %>% rename_with(~ paste0(.x, "_2"), -file),
    by = c("file_2" = "file"))

fwrite(dt_cor_count, sep='\t',
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/correlation_of_synonymous.tsv.gz")

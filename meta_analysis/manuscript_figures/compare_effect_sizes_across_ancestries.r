# Next, evaluate the Beta_Burden sumstats including and excluding Europe and compare the effect sizes

# First, let's get this working on the cluster - everything here up to reading the data back in is carried out on the cluster
# Read in all of the information

files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene", full.names=TRUE, recursive=TRUE)
gene_files <- grep("/n_cases_100/", files, value=TRUE)
gene_files <- gene_files[!grepl("extra_cauchy", gene_files)]

# First use the ALL of EUR file to define the set of pairs that make it through (this'll keep the file size down).
gene_files_all <- gene_files[!grepl(".*cutoff.([A-Za-z_]*).tsv.gz", gene_files)]

gene_files_all_list <- list()
for (file in gene_files_all)
{
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	cat(paste0(phe, "\n"))
	gene_files_all_list[[phe]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted", Pvalue < 6.7e-7) %>%
		mutate(phenotype = phe, ancestry = "ALL")
	setkeyv(gene_files_all_list[[phe]], c("Region", "Group", "max_MAF", "phenotype"))
}

gene_files_pop <- setdiff(gene_files, gene_files_all)
# Loop over the ancestries and loop over the phenotypes, merging at each step.
gene_files_list <- list()
for (file in gene_files_pop) {
	# Extract the phenotype name
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	ancestry <- ifelse(grepl(".*cutoff.([A-Za-z_]*).tsv.gz", file), gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file), "ALL")
	cat(paste0(phe, ": ", ancestry, "\n"))
	if (is.null(gene_files_list[[phe]])) {
		gene_files_list[[phe]] <- list()
	}
	gene_files_list[[phe]][[ancestry]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted") %>%
		mutate(phenotype = phe, ancestry = gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)) %>%
		rename(
			Pvalue_pop = Pvalue,
			BETA_Burden_pop = BETA_Burden,
			chisq_het_pop = chisq_het,
			Pvalue_het_pop = Pvalue_het,
			SE_Burden_pop = SE_Burden,
			BETA_meta_pop = BETA_meta,
			sum_weights_pop = sum_weights,
			ancestry_pop = ancestry,
			df_pop = df 
			) %>% select(-c("Stat", "type", "class")) 
	setkeyv(gene_files_list[[phe]][[ancestry]], c("Region", "Group", "max_MAF", "phenotype"))
	gene_files_list[[phe]][[ancestry]] <- merge(gene_files_list[[phe]][[ancestry]], gene_files_all_list[[phe]] %>% select(-c("Stat", "type", "class")))
}

for (phe in names(gene_files_list)) {
	gene_files_list[[phe]] <- rbindlist(gene_files_list[[phe]])
}
gene_files_list <- rbindlist(gene_files_list)
gene_files_list <- gene_files_list %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

fwrite(gene_files_list, file="meta_analysis/Burden_comparison_for_plotting_vs_ALL.tsv.gz", sep='\t')
gene_files_list <- fread("meta_analysis/Burden_comparison_for_plotting_vs_EUR.tsv.gz")
ggplot(gene_files_list, aes(x=BETA_Burden_pop, y=BETA_Burden, col=Group)) + geom_point() + facet_wrap(~ancestry_pop+max_MAF+case_control)

# The same thing, but comparing to EUR meta

# First use the ALL of EUR file to define the set of pairs that make it through (this'll keep the file size down).
gene_files_EUR <- gene_files[grepl(".*cutoff.EUR.tsv.gz", gene_files)]

gene_files_EUR_list <- list()
for (file in gene_files_EUR)
{
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	cat(paste0(phe, "\n"))
	gene_files_EUR_list[[phe]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted", Pvalue < 6.7e-7) %>%
		mutate(phenotype = phe, ancestry = "EUR")
	setkeyv(gene_files_EUR_list[[phe]], c("Region", "Group", "max_MAF", "phenotype"))
}

gene_files_pop <- setdiff(gene_files, c(gene_files_all, gene_files_EUR))
# Loop over the ancestries and loop over the phenotypes, merging at each step.
gene_files_list <- list()
for (file in gene_files_pop) {
	# Extract the phenotype name
	phe <- gsub(".*/([A-Za-z0-9]+)_.*", "\\1", file)
	ancestry <- gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)

	cat(paste0(phe, ": ", ancestry, "\n"))
	if (is.null(gene_files_list[[phe]])) {
		gene_files_list[[phe]] <- list()
	}

	gene_files_list[[phe]][[ancestry]] <- fread(file) %>% 
		filter(class == "Burden", type == "Inverse variance weighted") %>%
		mutate(phenotype = phe, ancestry = gsub(".*cutoff.([A-Za-z_]*).tsv.gz", "\\1", file)) %>%
		rename(
			Pvalue_pop = Pvalue,
			BETA_Burden_pop = BETA_Burden,
			chisq_het_pop = chisq_het,
			Pvalue_het_pop = Pvalue_het,
			SE_Burden_pop = SE_Burden,
			BETA_meta_pop = BETA_meta,
			sum_weights_pop = sum_weights,
			ancestry_pop = ancestry,
			df_pop = df 
			) %>% select(-c("Stat", "type", "class")) 
	setkeyv(gene_files_list[[phe]][[ancestry]], c("Region", "Group", "max_MAF", "phenotype"))
	gene_files_list[[phe]][[ancestry]] <- merge(gene_files_list[[phe]][[ancestry]], gene_files_EUR_list[[phe]] %>% select(-c("Stat", "type", "class")))
}

for (phe in names(gene_files_list)) {
	gene_files_list[[phe]] <- rbindlist(gene_files_list[[phe]])
}
gene_files_list <- rbindlist(gene_files_list)
gene_files_list <- gene_files_list %>% mutate(case_control =
	ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

fwrite(gene_files_list, file="data/meta_analysis/Burden_comparison_for_plotting_vs_EUR.tsv.gz", sep='\t')
gene_files_list <- fread("data/meta_analysis/Burden_comparison_for_plotting_vs_EUR.tsv.gz") %>% mutate(
	sign_BETA_Burden = sign(BETA_Burden),
	BETA_Burden = sign_BETA_Burden * BETA_Burden,
	BETA_Burden_pop = sign_BETA_Burden * BETA_Burden_pop)

# Factor to get the names correct
gene_files_list <- gene_files_list %>% 
	mutate(ancestry_pop = ifelse(ancestry_pop == "non_EUR", "non-EUR", ancestry_pop)) %>%
	mutate(ancestry_pop = factor(ancestry_pop, levels = c("AFR", "AMR", "EAS", "SAS", "non-EUR")))

pdf(width = 5.5, height = 3.5, file = "effect_sizes_binary.pdf")
ggplot(gene_files_list %>% filter(case_control, max_MAF != "0.01", Group == "pLoF", SE_Burden_pop < 0.06), aes(x=BETA_Burden_pop, y=BETA_Burden, col=Group)) + 
	geom_point() + facet_grid(max_MAF ~ ancestry_pop) + theme_minimal() + 
	scale_color_manual(values = colors) + geom_abline(intercept = 0 , slope = 1, colour = 'grey40', linetype = "dashed") +
	labs(x = TeX("$\\beta_{Burden}"), y = TeX("$\\beta_{Burden}$ (EUR)")) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf(width = 5.5, height = 3.5, file = "effect_sizes_cts.pdf")
ggplot(gene_files_list %>% filter(!case_control, max_MAF != "0.01", Group == "pLoF", SE_Burden_pop < 0.06), aes(x=BETA_Burden_pop, y=BETA_Burden, col=Group)) + 
	geom_point() + facet_grid(max_MAF~ancestry_pop) + theme_minimal() + 
	scale_color_manual(values = colors) + geom_abline(intercept = 0 , slope = 1, colour = 'grey40', linetype = "dashed") +
	labs(x = TeX("$\\beta_{Burden}"), y = TeX("$\\beta_{Burden}$ (EUR)")) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(width = 7.5, height = 3, file = "effect_sizes_binary.pdf")
ggplot(gene_files_list %>% filter(case_control, max_MAF == "0.001", Group == "pLoF", SE_Burden_pop < 0.06), aes(x=BETA_Burden_pop, y=BETA_Burden, col=Group)) + 
	geom_point() + facet_grid(~ancestry_pop) + theme_minimal() + 
	scale_color_manual(values = colors) + geom_abline(intercept = 0 , slope = 1, colour = 'grey40', linetype = "dashed") +
	labs(x = TeX("$\\beta_{Burden}"), y = TeX("$\\beta_{Burden}$ (EUR)")) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
pdf(width = 7.5, height = 3, file = "effect_sizes_cts.pdf")
ggplot(gene_files_list %>% filter(!case_control, max_MAF == "0.001", Group == "pLoF", SE_Burden_pop < 0.06), aes(x=BETA_Burden_pop, y=BETA_Burden, col=Group)) + 
	geom_point() + facet_grid(~ancestry_pop) + theme_minimal() + 
	scale_color_manual(values = colors) + geom_abline(intercept = 0 , slope = 1, colour = 'grey40', linetype = "dashed") +
	labs(x = TeX("$\\beta_{Burden}"), y = TeX("$\\beta_{Burden}$ (EUR)")) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
	theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

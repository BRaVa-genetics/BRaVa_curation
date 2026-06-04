library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)

source("../meta_analysis_utils.r")

gene_files_list <- fread(paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
		"Burden_comparison_for_plotting_vs_EUR.tsv.gz")) %>% 
	mutate(
		sign_BETA_Burden = sign(BETA_Burden),
		BETA_Burden = sign_BETA_Burden * BETA_Burden,
		BETA_Burden_pop = sign_BETA_Burden * BETA_Burden_pop,
		max_MAF = as.numeric(max_MAF))

dt_deming <- fread("../manuscript_tables/Tables/deming_regressions.tsv.gz") %>%
	rename(ancestry_pop = ancestry,
		max_MAF = `max MAF`)

# Factor to get the names correct
gene_files_list <- gene_files_list %>% 
	mutate(ancestry_pop = ifelse(
			ancestry_pop == "non_EUR", "non-EUR", ancestry_pop)) %>%
	mutate(ancestry_pop = factor(ancestry_pop,
			levels = c("AFR", "AMR", "EAS", "SAS", "non-EUR"))) %>%
	mutate(Group = case_when(
		Group == "damaging_missense_or_protein_altering" ~ "DM/PA",
		Group == "pLoF;damaging_missense_or_protein_altering" ~ "pLoF;DM/PA",
		TRUE ~ Group))

pdf(width = 5.5, height = 3.5, file = "Figures/effect_size_comparisons.pdf")
for (cc in c(TRUE, FALSE)) {
	for (g in c("pLoF", "DM/PA", "pLoF;DM/PA")) {
		p <- ggplot(gene_files_list %>% 
			filter(case_control == !!cc,
				max_MAF != 0.01,
				Group == !!g,
				SE_Burden_pop < 0.06),
			aes(x=BETA_Burden_pop, y=BETA_Burden, col=ancestry_pop)) + 
		geom_point() + facet_grid(max_MAF ~ ancestry_pop) + theme_minimal() + 
		scale_color_manual(values=pop_colors) + 
		geom_abline(intercept = 0 , slope = 1,
			colour = 'grey40', linetype = "dashed") +
		labs(x = TeX("$\\beta_{Burden}"),
			y = TeX("$\\beta_{Burden}$ (EUR)")) + 
		geom_vline(xintercept=0) + geom_hline(yintercept=0) +
		theme(legend.position = "none") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1))

		p <- p + geom_abline(
			data = dt_deming %>% filter(
				case_control == !!cc,
				Group == !!g),
			aes(intercept = Intercept, slope = Slope),
			color = "black"
		)
		print(p)
	}
}
dev.off()

pdf(width = 7.5, height = 3, file = "Figures/effect_sizes_0.001.pdf")
for (cc in c(TRUE, FALSE)) {
	for (g in  c("pLoF", "DM/PA", "pLoF;DM/PA")) {
		p <- ggplot(gene_files_list %>%
			filter(case_control == !!cc,
				max_MAF == 0.001,
				Group == !!g,
				SE_Burden_pop < 0.06),
			aes(x=BETA_Burden_pop, y=BETA_Burden, col=ancestry_pop)) + 
		geom_point() + facet_grid(~ancestry_pop) + theme_minimal() + 
		scale_color_manual(values = pop_colors) +
		geom_abline(intercept = 0 , slope = 1,
			colour = 'grey40', linetype = "dashed") +
		labs(x = TeX("$\\beta_{Burden}"),
			y = TeX("$\\beta_{Burden}$ (EUR)")) + 
		geom_vline(xintercept=0) + geom_hline(yintercept=0) +
		theme(legend.position = "none") +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))

		dt_sub <- dt_deming %>% filter(
				case_control == !!cc,
				Group == !!g,
				max_MAF == 0.001)
		print(cc)
		print(g)
		print(dt_sub)

		if (nrow(dt_sub) > 0) {
			p <- p + geom_abline(
				data = dt_sub,
				aes(intercept = Intercept, slope = Slope),
				color = "black"
			)
		}
		print(p)
	}
}
dev.off()

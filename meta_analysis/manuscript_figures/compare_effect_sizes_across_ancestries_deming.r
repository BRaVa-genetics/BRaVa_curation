# We have plotting comparing effect sizes but not an explicit test
# Here we create those collections of tests

# Take a look at Jurgens to see what they did

# Deming regression

# Get this sorted. They definitely included a bunch of caveats to get it to work

library(deming)
library(dplyr)
library(data.table)

gene_files_list <- fread("Burden_comparison_for_plotting_vs_EUR.tsv.gz") %>% 
	mutate(
		sign_BETA_Burden = 1, #sign(BETA_Burden),
		BETA_Burden = sign_BETA_Burden * BETA_Burden,
		BETA_Burden_pop = sign_BETA_Burden * BETA_Burden_pop)

gene_files_list <- gene_files_list %>% 
	mutate(ancestry_pop = ifelse(
			ancestry_pop == "non_EUR", "non-EUR", ancestry_pop)) %>%
	mutate(ancestry_pop = factor(ancestry_pop,
			levels = c("AFR", "AMR", "EAS", "SAS", "non-EUR"))) %>%
	mutate(Group = case_when(
		Group == "damaging_missense_or_protein_altering" ~ "DM/PA",
		Group == "pLoF;damaging_missense_or_protein_altering" ~ "pLoF;DM/PA",
		TRUE ~ Group))

dt_deming <- data.table(
	case_control = logical(),
	ancestry = character(),
	`max MAF` = numeric(),
	Group = character(),
	Slope = numeric(),
	`Slope SE` = numeric(),
	`Slope P-value` = numeric(),
	# `Slope P-value != 1` = numeric(),
	`Slope 95% CI` = character(), 
	Intercept = numeric(),
	`Intercept SE` = numeric(),
	`Intercept P-value` = numeric(),
	`N pairs` = integer())

for (cc in c(TRUE, FALSE)) {
	for (g in  c("pLoF", "DM/PA", "pLoF;DM/PA")) {
		for (anc in c("AFR", "AMR", "EAS", "SAS", "non-EUR")) {
			for (mm in c(0.001, 0.0001)) {

				dt <- gene_files_list %>% filter(
					ancestry_pop == !!anc,
					case_control == !!cc,
					Group == !!g,
					max_MAF == !!mm,
					!is.na(BETA_meta),
					!is.na(BETA_meta_pop))

				fit <- deming(BETA_meta_pop ~ BETA_meta, 
					xstd = SE_Burden,
					ystd = SE_Burden_pop,
					data = dt)

				slope <- fit$coef[2]
				se_slope <- fit$se[2]
				z_stat <- slope / se_slope
				log10_pval <- log10(2) + pnorm(-abs(z_stat), log.p=TRUE)/log(10)
				cat(anc, g, ifelse(cc, "case-control", "cts"), mm, "\n")
				cat("slope:", slope, "\n")
				cat("se slope", se_slope, "\n")
				cat("log10(p-value):", log10_pval, "\n")
				cat("number of observations", nrow(dt), "\n")

				# z_stat_1 = (slope - 1) / se_slope
				# log10_pval_1 <- log10(2) + pnorm(-abs(z_stat_1), log.p=TRUE)/log(10)
				dt_deming <- rbind(
					dt_deming,
					data.table(
						case_control = cc,
						ancestry = anc,
						`max MAF` = mm,
						Group = g,
						Slope = slope,
						`Slope SE` = se_slope,
						`Slope P-value` = log10_pval,
						# `Slope P-value != 1` = log10_pval_1,
						`Slope 95% CI` = paste0("[",
							format(round(slope - 1.96*se_slope, 2), nsmall = 2), ", ",
							format(round(slope + 1.96*se_slope, 2), nsmall = 2), "]"),
						Intercept = fit$coef[1],
						`Intercept SE` = fit$se[1],
						`Intercept P-value` = (
							pnorm(
								-abs(fit$coef[1]/fit$se[1]),
								log.p=TRUE)/log(10)),
						`N pairs` = nrow(dt)
					)
				)
			}
		}
	}
}

dt_deming <- dt_deming %>% filter(`N pairs` > 40)

fwrite(dt_deming, sep='\t', quote=FALSE, file="../manuscript_tables/Tables/deming_regressions.tsv.gz")

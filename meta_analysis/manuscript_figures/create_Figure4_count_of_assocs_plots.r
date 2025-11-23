#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2)

source("../meta_analysis_utils.r")

# Copy this file down from BMRC to create plots.
# BRMC location: /well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_data.tsv.gz
# Ensure the above file is updated following rounds of meta-analysis and incorporation of additional data.
# scp qen698@cluster2.bmrc.ox.ac.uk:/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_data.tsv.gz ~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/

plot_unique <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_data.tsv.gz")
# Determine the counts of binary and continuous traits in the results:
files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100", full.names=TRUE)
files <- grep("cutoff.tsv.gz", files, value=TRUE)
phenos <- unique(gsub(".*/([A-Za-z0-9]+)_.*", "\\1", files))
n_binary <- sum(phenos %in% phenotype_class$binary)
n_cts <- sum(phenos %in% phenotype_class$continuous)
# Ensure that sex specific traits and 'ALL' traits are restricted to those,
# otherwise CCPM will have multiple shots at hitting a significant association

for (cc in c(TRUE, FALSE)) {
	pdf(width=ifelse(cc, 5, 3.5), height=4, file=ifelse(cc, "Figures/binary_unique_hits_count.pdf", "Figures/cts_unique_hits_count.pdf"))
	plot_unique <- plot_unique %>%
  		mutate(dataset = factor(dataset, levels = c(setdiff(sort(unique(dataset)), "Meta"), "Meta")))
	p <- ggplot(plot_unique %>% filter(case_control == !!cc),
				aes(x=dataset, y=count, fill=ancestry)) + 
			geom_bar(stat= "identity", position = "dodge") + theme_minimal() +
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
			scale_y_continuous(labels = scales::comma) +
			scale_fill_manual(values=pop_colors) + 
			labs(x = NULL, y = "Number of unique\n(gene, trait) associations", fill = "Genetic\nancestry",
				title=ifelse(cc,
					paste0("Binary traits (N=", n_binary, ")"),
					paste0("Continuous traits (N=", n_cts, ")")
					)
				)
	print(p)
	dev.off()
}

plot_unique <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/plot_unique_hits_no_height_data.tsv.gz")
n_cts <- n_cts - 1
# Ensure that sex specific traits and 'ALL' traits are restricted to those,
# otherwise CCPM will have multiple shots at hitting a significant association

cc <- FALSE
pdf(width=ifelse(cc, 5, 3.5), height=4, file="Figures/cts_unique_hits_no_height_count.pdf"))
plot_unique <- plot_unique %>%
  	mutate(dataset = factor(dataset, levels = c(setdiff(sort(unique(dataset)), "Meta"), "Meta")))
p <- ggplot(plot_unique %>% filter(case_control == !!cc),
		aes(x=dataset, y=count, fill=ancestry)) + 
		geom_bar(stat= "identity", position = "dodge") + theme_minimal() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
		scale_y_continuous(labels = scales::comma) +
		scale_fill_manual(values=pop_colors) + 
		labs(x = NULL, y = "Number of unique\n(gene, trait) associations", fill = "Genetic\nancestry",
			title=ifelse(cc,
				paste0("Binary traits (N=", n_binary, ")"),
				paste0("Continuous traits (N=", n_cts, ")")
			)
		)
print(p)
dev.off()

library(data.table)
library(dplyr)
library(ggplot2)

source("../meta_analysis_utils.r")

meta_files <- grep(".vcf.gz$",
	dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100",
		full.names=TRUE), value=TRUE)
dt_plot_list <- list()
for (file in meta_files)
{
	cat(file, "\n")
	cmd <- paste(
		"bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", file)
	dt <- fread(cmd = cmd) %>% 
		rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7,
			`P-value`=V8) %>%
		mutate(
			BETA=as.numeric(BETA), SE=as.numeric(SE),
			`P-value`=-as.numeric(`P-value`))
	dt <- data.table(dt) %>% filter(!is.na(`P-value`))
	setkey(dt, "P-value")
	dt_plot_list[[file]] <- dt
	phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", file)
	type <- ifelse(phenotype %in% phenotype_class$continuous, "continuous", "binary")
	dt_plot_list[[file]]$phenotype <- phenotype
	dt_plot_list[[file]]$type <- type
}

dt_manhattan_plots <- rbindlist(dt_plot_list)
fwrite(dt_manhattan_plots,
	file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/manhattan_plots.tsv.gz",
	sep="\t", quote=FALSE)

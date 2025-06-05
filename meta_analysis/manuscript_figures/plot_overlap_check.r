library(data.table)
library(dplyr)
library(ggplot2)
library(hexbin)
library(latex2exp)

source("meta_analysis_utils.r")

sample_overlap_folder <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100_sample_overlap"
no_sample_overlap_folder <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100_no_overlap"


files <- grep("gz$", dir(sample_overlap_folder), value=TRUE)
files_sample_overlap <- paste0(sample_overlap_folder, "/", files)
files_no_sample_overlap <- paste0(no_sample_overlap_folder, "/", files)

pdf("overlap_check_binned.pdf", width=3, height=3)
for (i in 1:length(files_sample_overlap)) {
    print(i)
	dt2 <- fread(files_sample_overlap[i])
	dt1 <- fread(files_no_sample_overlap[i])
	setkeyv(dt1, c("Region", "Group", "max_MAF", "class", "type"))
	setkeyv(dt2, c("Region", "Group", "max_MAF", "class", "type"))
	dt <- merge(dt1, dt2)
	phenotype <- renaming_phenotype_list[gsub(".*/([A-Za-z0-9]+)_.*", "\\1", files_sample_overlap[i])]
	p <- ggplot(dt, aes(x=-log10(`Pvalue.x`), y=-log10(`Pvalue.y`))) + geom_hex(bins=100) +
	labs(x=TeX("$-\\log_{10}(P_{independent})$"), y=TeX("$-\\log_{10}(P_{overlap})$"),
		title = phenotype) + theme_minimal() +
	geom_abline(slope = 1, intercept = 0, color = "red") + theme(legend.position = "none")
	print(p)
}
dev.off()

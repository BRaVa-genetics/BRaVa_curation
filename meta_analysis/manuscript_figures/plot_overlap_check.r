library(data.table)
library(dplyr)
library(ggplot2)
library(hexbin)
library(latex2exp)

source("../meta_analysis_utils.r")

sample_overlap_folder <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100_sample_overlap"
no_sample_overlap_folder <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100"

files <- grep("gz$", dir(sample_overlap_folder), value=TRUE)
files_sample_overlap <- paste0(sample_overlap_folder, "/", files)
files_no_sample_overlap <- paste0(no_sample_overlap_folder, "/", files)
n_bins <- 100

pdf("Figures/overlap_check_binned.pdf", width=4, height=3)
dt_list <- list()
for (i in 1:length(files_sample_overlap)) {
    print(i)
	dt2 <- fread(files_sample_overlap[i]) %>% filter(Group != "Cauchy", max_MAF != "0.01")
	dt1 <- fread(files_no_sample_overlap[i]) %>% filter(Group != "Cauchy", max_MAF != "0.01")
	setkeyv(dt1, c("Region", "Group", "max_MAF", "class", "type"))
	setkeyv(dt2, c("Region", "Group", "max_MAF", "class", "type"))
	dt_list[[i]] <- merge(dt1, dt2)
	eps <- 1e-300
	dt_list[[i]][, Pvalue_overlap := -log10(pmax(Pvalue.y, eps))]
	dt_list[[i]][, `Pvalue` := -log10(pmax(`Pvalue.x`, eps))]
	phenotype <- renaming_phenotype_list[gsub(".*/([A-Za-z0-9]+)_.*", "\\1", files_sample_overlap[i])]
	p <- ggplot(dt_list[[i]], aes(x = Pvalue, y = Pvalue_overlap)) + geom_hex(bins=n_bins) +
		labs(x=TeX("$-\\log_{10}(P)$"), y=TeX("$-\\log_{10}(P_{overlap})$"),
			title = phenotype) + theme_minimal(base_size = 12) +
	  	scale_fill_viridis_c(trans = "log10", guide = guide_colorbar(title = "")) +
	  	theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11)) +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", col='indianred3')
	print(p)
}
dev.off()

dt <- rbindlist(dt_list, fill=TRUE)

pdf("Figures/overlap_check_binned_combined.pdf", width=4, height=3)
p <- ggplot(dt, aes(x = Pvalue, y = Pvalue_overlap)) + geom_hex(bins=n_bins) +
	labs(x=TeX("$-\\log_{10}(P)$"), y=TeX("$-\\log_{10}(P_{overlap})$")) + 
	theme_minimal(base_size = 12) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", col='indianred3') + 
	scale_fill_viridis_c(trans = "log10", guide = guide_colorbar(title = "")) +
	theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11))
print(p)

dt_zoom <- dt[Pvalue <= 10 & Pvalue_overlap <= 10]
p <- ggplot(dt_zoom, aes(x = Pvalue, y = Pvalue_overlap)) + geom_hex(bins=n_bins) +
	labs(x=TeX("$-\\log_{10}(P)$"), y=TeX("$-\\log_{10}(P_{overlap})$")) + 
	theme_minimal(base_size = 12) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed", col='indianred3') + 
	scale_fill_viridis_c(trans = "log10", guide = guide_colorbar(title = "")) +
	theme(panel.grid.minor = element_blank(), strip.text = element_text(size = 11))
print(p)
dev.off()

dt_signif <- dt %>% filter((Pvalue > -log10(2.5e-7) & class != "Burden") | (Pvalue > -log10(6.7e-7) & class == "Burden"))
dt_signif <- dt_signif %>% mutate(retained = (Pvalue_overlap > -log10(2.5e-7) & class != "Burden") | (Pvalue_overlap > -log10(6.7e-7) & class == "Burden"))

dt %>% filter(Group %in% c("pLoF", "damaging_missense_or_protein_altering", "pLoF;damaging_missense_or_protein_altering"))

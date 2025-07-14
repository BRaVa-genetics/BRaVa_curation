library(data.table)
library(dplyr)
library(ggplot2)

length <- 10000
maxP <- log10(1)
ribbon_p <- 0.95
source("../QC/utils/pretty_plotting.r")
source("../phenotypes/BRaVa_phenotypes_utils.r")
source("meta_analysis_utils.r")

pdf(file="manuscript_figures/Figures/QQ_variant.pdf", width=6, height=4)
meta_files <- grep(".vcf.gz$", dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100", full.names=TRUE), value=TRUE)

for (file in meta_files)
{
	cat(file, "\n")
	cmd <- paste("bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", file)
	dt <- fread(cmd = cmd) %>% 
		rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7,
			`P-value`=V8) %>%
		mutate(BETA=as.numeric(BETA), SE=as.numeric(SE),
			`P-value`=-as.numeric(`P-value`))
	dt <- data.table(dt) %>% filter(!is.na(`P-value`))
	setkey(dt, "P-value")
	phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", file)
	phenotype <- renaming_phenotype_list[[phenotype]]
	
	# We want the QQ to be linear on the log10 scale, so we should sample in 
	# that way
	nP <- nrow(dt)
	minP <- log10(1/(nP+1))
	elements <- unique(round(10^seq(minP, maxP, length.out=length) * nP))
	# Ensure that elements contains the top 100 associations
	elements <- sort(union(elements, seq(1,100)))

	dt_plot <- dt[elements, ]
	dt_plot[ , P_expected := -log10(seq(10^minP, 10^maxP,
		length.out=(nP+1))[elements])]
	dt_plot <- dt_plot %>% 
		rename(P_observed = `P-value`) %>% 
		mutate(
			P_observed = ifelse(P_observed < -320, 320, -P_observed),
			clower = -log10(qbeta(p = (1 - ribbon_p) / 2,
				shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			cupper = -log10(qbeta(p = (1 + ribbon_p) / 2,
				shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			OR = exp(BETA)
		)

	type <- ifelse(phenotype %in% phenotype_class$continuous,
		"continuous", "binary")

	if (type == "continuous") {
		dt_plot$color <- cut(dt_plot$BETA,
        breaks = c(-Inf, -0.5, 0, 0.5, Inf),
        labels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    dt_plot$color <- factor(dt_plot$color,
	        levels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    dummy_data <- data.frame(P_expected = NA, P_observed = NA,
	        color = factor(c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"),
	        levels = c("< -0.5", "[-0.5, 0)", "[0, 0.5]", "> 0.5"))
	    )
	} else {
		dt_plot$color <- cut(dt_plot$OR,
        breaks = c(-Inf, 0.5, 1, 2, Inf),
        labels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    dt_plot$color <- factor(dt_plot$color,
	        levels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    dummy_data <- data.frame(P_expected = NA, P_observed = NA,
	        color = factor(c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"),
	        levels = c("< 0.5", "[0.5, 1)", "[1, 2]", "> 2"))
	    )
	}

	p <- create_pretty_qq_plot(
		plot_title=phenotype,
		plot_subtitle="coding regions; gnomAD popmax < 0.01",
		rbind(dt_plot, dummy_data, fill=TRUE),
		aes(x=P_expected, y=P_observed, col=color), key_cols=c("P_observed"),
		aes_ribbon = aes(ymin=clower, ymax=cupper),
		x_label=TeX("$-\\log_{10}(P_{expected})$"), 
		y_label=TeX("$-\\log_{10}(P_{observed})$"),
		print_p=FALSE,
		)

	if (type == "continuous")
	{
		p <- p + scale_color_manual(
            values = c(
                "< -0.5" = "blue3",
                "[-0.5, 0)" = "cornflowerblue",
                "[0, 0.5]" = "indianred3",
                "> 0.5" = "red"),
            labels = c("< -0.5" = "< -0.5",
                "[-0.5, 0)" = "[-0.5, 0)",
                "[0, 0.5]" = "[0, 0.5]",
                "> 0.5" = "> 0.5"),
            name = "Effect size",  aesthetics = c("colour", "fill")
            ) + guides(colour = guide_legend(override.aes = list(size=5)))
	} else {
		p <- p + scale_color_manual(
			values = c(
				"< 0.5" = "blue3",
                "[0.5, 1)" = "cornflowerblue",
                "[1, 2]" = "indianred3",
                "> 2" = "red"),
            labels = c("< 0.5" = "< 0.5",
                "[0.5, 1)" = "[0.5, 1)",
                "[1, 2]" = "[1, 2]",
                "> 2" = "> 2"),
            name = "Odds ratio",  aesthetics = c("colour", "fill")
            ) + guides(colour = guide_legend(override.aes = list(size=5)))
	}
	print(p)
}
dev.off()

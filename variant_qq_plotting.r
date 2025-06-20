library(data.table)
library(dplyr)
library(ggplot2)

length <- 10000
maxP <- log10(1)
ribbon_p <- 0.95
source("../QC/utils/pretty_plotting.r")
source("../phenotypes/BRaVa_phenotypes_utils.r")

phenotype_class <- list(
	binary = c("AAA","AcApp","AcuLymLeuk","Adenomy","AMD","ALamy","AUD","AloAre",
		"AnoNer","AoSten","Asth","AtopDis","AFib","ADHD","ASD","BCLL","BenCervUterNeo",
		"BenIntNeo","BenNodGoit","BladCanc","BrainCNSCanc","BreastCanc","BrugSynd",
		"BuliNer","BullPemph","CarShock","HCM","CRVO","CervCanc","CML","COPD","CRF",
		"CoffSirSynd","ColonRectCanc","CAD","CCANS","EatDis","Endocar","Endometr",
		"EsophCanc","EssThrom","EFRMB","FSP","FemInf","FemInfAC","FolLymph","Gout",
		"GravesDis","HemoChromo","HF","HepCarcin","HTN","HHD","HypoThyr","HypoThyrSec",
		"IPF","ITP","IBD","IFHern","ILDSarc","IodDef","KabSynd","KidCanc","KleefSynd",
		"LaryxCanc","Leuk","LiverCanc","LiverFibCirr","LongQTSynd","LymphThyrit",
		"MalInf","MatHem","MatHypDis","MODYDiab","MultiMyel","MS","MECS","Myocard",
		"Narco1","NonFuncPitAd","NHL","NonPapTCCBlad","NonRheuValv","OUD","OCD",
		"OvCanc","Pancreat","ParkDis","PeptUlcer","PAD","PlacInsuf","PCOS",
		"PolycythVera","Preeclamps","PregLoss","POAG","PrimSjoSynd","Prolactinom",
		"Psori","RheumHeaDis","RheumArth","RomWardSynd","Sarcoid","SebDerm",
		"SpinaBifAp","StomCanc","Stroke","SLE","TAAD","ThyroCanc","T2Diab","Urolith",
		"UterCanc","VaricVeins","VTE"),
	continuous = c("ALT", "AlcCons", "AST","BMI","CRP","CACS","CK","HDLC","Height",
		"LDLC","TChol","TG","WHRBMI", "LVH","Append","HipRep","CogAbil","EduAtt",
		"PsySymp","SchGrades","SCDCAT")
	)

for (phe in BRaVa_pilot_phenotypes)
{
	cat(phe, "\n")
	files <- grep(paste0(phe, ".*vcf.gz$"),
			dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/vcf",
				recursive=TRUE, full.names=TRUE), value=TRUE)
	meta_files <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/n_cases_100/", full.names)
	meta_files <- grep("vcf.gz$", meta_files, value=TRUE)

	meta_file <- grep(paste0(phe, "_"), meta_files, value=TRUE)
	if (length(meta_file) == 1) {
		cmd <- paste("bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", meta_file)
		dt <- fread(cmd = cmd) %>% 
			rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7, `P-value`=V8) %>%
			mutate(BETA=as.numeric(BETA), SE=as.numeric(SE), `P-value`=-as.numeric(`P-value`))
		dt_meta <- dt %>% filter(is.na(`P-value`))
		setkeyv(dt_meta, c("CHR", "POS", "REF", "ALT"))
		if (nrow(dt_meta) > 0) {
			for (file in files) {
				print(file)
				cmd <- paste("bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", file)
				dt <- fread(cmd = cmd) %>% 
					rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7, `P-value`=V8) %>%
					mutate(BETA=as.numeric(BETA), SE=as.numeric(SE), `P-value`=-as.numeric(`P-value`))
				setkeyv(dt, c("CHR", "POS", "REF", "ALT"))
				dt_merge <- merge(dt, dt_meta)
				if (any(dt_merge$`BETA.x` == 0)) {
					print(dt_merge)
				}
			}
		}
	} else {
		cat("meta-analysis results file for ", phe, "doesn't exist...what's going on?\n")
		print(meta_file)
	}
}


pdf(file="QQ_variant.pdf", width=6, height=4)
for (file in files)
{
	# dt <- fread(file, key="P-value")
	cmd <- paste("bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\t[ %ES]\t[ %SE]\t[ %LP]\n'", file)
	dt <- fread(cmd = cmd) %>% 
		rename(ID=V1, CHR=V2, POS=V3, REF=V4, ALT=V5, BETA=V6, SE=V7, `P-value`=V8) %>%
		mutate(BETA=as.numeric(BETA), SE=as.numeric(SE), `P-value`=-as.numeric(`P-value`))
	dt <- data.table(dt) %>% filter(!is.na(`P-value`))
	setkey(dt, "P-value")
	phenotype <- gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", file)
	
	# We want the QQ to be linear on the log10 scale, so we should sample in that way
	nP <- nrow(dt)
	minP <- log10(1/(nP+1))
	elements <- unique(round(10^seq(minP, maxP, length.out=length) * (nP+1)))
	# Ensure that elements contains the top 100 associations
	elements <- sort(union(elements, seq(1,100)))

	dt_plot <- dt[elements, ]
	dt_plot[ , P_expected := -log10(seq(10^minP, 10^maxP, length.out=(nP+1))[elements])]
	dt_plot <- dt_plot %>% 
		rename(P_observed = `P-value`) %>% 
		mutate(
			P_observed = ifelse(P_observed < -320, 320, -P_observed),
			clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = (nP:1)[elements], shape1 = (1:nP)[elements])),
			OR = exp(BETA)
		)

	type <- ifelse(phenotype %in% phenotype_class$continuous, "continuous", "binary")

	if (type == "continuous") {
		dt_plot$color <- cut(dt_plot$Effect,
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

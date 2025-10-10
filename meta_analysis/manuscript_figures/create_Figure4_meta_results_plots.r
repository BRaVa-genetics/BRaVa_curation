#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
devtools::install_github("mkanai/rgsutil")
library(rgsutil)
source("../meta_analysis_utils.r")
source("../../QC/utils/pretty_plotting.r")

significance_T1 <- 6.7e-7
significance_T2 <- 2.5e-7
loglog <- TRUE

get_category_colors <- function(category=NULL)
{
	color_infectious <- "#822B20"
	color_neoplasms <- "#D26529"
	color_immune <- "#C6B844"
	color_metabolic <- "#999999"
	color_behavioral <- "#B17D89"
	color_nervous <- "#77216F"
	color_eye <- "#97C1A9"
	color_ear <- "#1F2D58"
	color_circulatory <- "#4F572E"
	color_respiratory <- "#B87726"
	color_digestive <- "#83BAE4"
	color_dermatologic <- "#1657C4"
	color_musculoskeletal <- "#5891E0"
	color_urinary <- "#F2C4C2"
	color_pregnancy <- "#F3C44E"
	color_other_factors <- "#4D7F36"

	category_colors <- list(
		Cardiovascular = color_circulatory,
		`Sense organs`= color_eye,
		Respiratory = color_respiratory,
		Neoplasms = color_neoplasms,
		Genitourinary = color_urinary,
		Muscloskeletal = color_musculoskeletal,
		Gastrointestinal = color_digestive,
		Dermatological = color_dermatologic,
		`Endocrine/Metabolic` = color_metabolic,
		Pregnancy = color_pregnancy
		)

	if (!is.null(category)) {
		return(category_colors[category])
	}
	return(category_colors)
}

make_gene_manhattan_category_plot <- function(dt, buffer=100000000,
	chr_lengths=chr_lengths_38, significance_T1=6.7e-7, significance_T2=2.5e-7, ggplot_theme=theme_bw,
	save_figure=FALSE, file='file_out', scaling=1, width=230, height=100, 
	print_p=FALSE, loglog=TRUE)
{  	
    contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    start_end <- get_start_and_end(chr_lengths)
    dt_contigs <- data.frame(contig=contigs_, start=start_end$start,
    	end=start_end$end)

    # Now, we also want to shift by phenotype category
    if (!(all(c("Pvalue", "phenotype_category", "chromosome_name") %in% names(dt)))) {
    	cat("Pvalue, phenotype_category, or chromosome_name not present in the passed data.table\n")
    }

    if (!("position" %in% names(dt))) {
    	dt <- dt %>% mutate(position = (
    		start_position + (end_position - start_position)/2))
    }

    # Take the vector of phenotype categories
    # Order them alphabetically
    phenotype_categories <- sort(unique(dt$phenotype_category))
    print(phenotype_categories)
    phenotype_category_list <- list()
    i <- 1
   	for (phenotype_category in phenotype_categories) {
   		phenotype_category_list[[phenotype_category]] <- i
   		i <- i+1
   	}
    dt_categories <- data.frame(
    	phenotype_category = phenotype_categories,
    	start = ((dt_contigs$end[nrow(dt_contigs)]) * (seq(1,length(phenotype_categories))-1) + 1),
    	end = ((dt_contigs$end[nrow(dt_contigs)]) * (seq(1,length(phenotype_categories))))) %>%
    	mutate(middle = floor(start + (end-start)/2), length = (end-start)) %>% 
    	mutate(shifted_position = middle + (unlist(phenotype_category_list[phenotype_category])-1) * buffer)

    dt_plot <- dt %>% transmute(
    	contig=gsub("chr", "", chromosome_name),
    	position=as.integer(position),
    	pval=as.numeric(Pvalue),
    	labels=external_gene_name,
    	phenotype_category=phenotype_category) %>% mutate(
    	x = dt_categories[unlist(phenotype_category_list[phenotype_category]), 'start'] + 
    		dt_contigs[gsub('X', 23, contig), 'start'] + 
    		position + (unlist(phenotype_category_list[phenotype_category])-1) * buffer)

  	dt_plot <- dt_plot %>% mutate(y = ifelse(pval < 1e-300, 300, -log10(pval)))

  	if (loglog) {
		transform_y <- function(y) {
			out <- numeric(length(y))
			log_part <- y <= 10
			out[log_part] <- y[log_part]
			out[!log_part] <- 10 + 2*(log2(y[!log_part] - 10 + 1))
			return(out)
		}
	} else {
		transform_y <- function(y) { return(y) }
	}

	# Apply transformation
	dt_plot <- dt_plot %>% mutate(y_trans = transform_y(y))
	breaks <- c(2, 4, 6, 8, 10, 20, 50, 100, 200, 300)
	breaks_trans <- transform_y(breaks)

    p <- ggplot(dt_plot, aes(x=x, y=y_trans, col=phenotype_category)) + geom_point_rast(size=0.5)
    if (loglog) {
	    p <- p + scale_y_continuous(
	    	breaks = breaks_trans,
	    	labels = breaks
	    )
    } else {
    	p <- p + scale_y_continuous(breaks=scales::pretty_breaks(n=10))
    }
    p <- p + geom_hline(yintercept=-log10(significance_T1), color='#E15759', linetype='dashed') +
    	geom_hline(yintercept=-log10(significance_T2), color='black', linetype='dashed') +
        scale_x_continuous(breaks=dt_categories$shifted_position, labels=dt_categories$phenotype_category) +
        scale_color_manual(values = get_category_colors()) +
        labs(x='', y=TeX("$-\\log_{10}(P)$")) + ggplot_theme() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    }

    if (print_p) { print(p) }

    return(list(p=p, dt=dt_plot))
}

if (loglog) {
	transform_y <- function(y) {
		out <- numeric(length(y))
		log_part <- y <= 10
		out[log_part] <- y[log_part]
		out[!log_part] <- 10 + 2*(log2(y[!log_part] - 10 + 1))
		return(out)
	}
} else {
	transform_y <- function(y) { return(y) }
}

file_root <- c("meta_analysis", "AFR", "AMR", "EAS", "EUR", "SAS", "non_EUR")

for (i in 1:length(file_root)) {
	meta_list <- fread(
		paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
			file_root[i], "_figure_4.tsv.gz"))
	meta_list <- meta_list %>% mutate(phenotype_category = unlist(phenotype_broad_categories[phenotype]))
	meta_list <- meta_list %>% filter(phenotype != "AlcCons")
	# Create the Burden, SKAT, and SKAT-O versions
	# Do the same thing, but split by case-control vs cts (way more power for cts).
	for (cc in c(TRUE, FALSE)) {
		cat(ifelse(cc, "case control\n", "continuous\n"))
		meta_list_tmp <- meta_list %>% filter(case_control == cc)
		p <- make_manhattan_plot(meta_list_tmp$chromosome_name,
			meta_list_tmp$start_position,
			meta_list_tmp$Pvalue,
			threshold=1000, significance_T = significance_T1,
			label=meta_list_tmp$external_gene_name, 
			colour_1 = "#6583E6",
			colour_2 = "#384980",
			loglog=TRUE)

		p$p <- p$p + geom_hline(yintercept=-log10(significance_T2), color='black', linetype='dashed')
		threshold <- ifelse(cc, 10, 50)

		if (nrow(p$dt %>% filter(y > transform_y(threshold))) < 20) {
			threshold <- 10
		}

		p$p <- p$p + geom_label_repel(
			data = unique(subset(p$dt, y > transform_y(threshold)) %>% group_by(labels) %>% 
				filter(y == max(y))) %>% ungroup(),
			size = 3, aes(label=labels),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',
			min.segment.length=0)
		width <- 230
		height <- 100
		scaling <- 1
		file <- paste0(file_root[i], "_",
			ifelse(cc, "unique_case_control", "unique_cts"))
		ggsave(paste0("Figures/", file, '.pdf'), p$p, width=width*scaling,
			height=height*scaling, units='mm')
		width <- 150
		p <- make_gene_manhattan_category_plot(meta_list_tmp, buffer=1000000000,
			scaling=scaling, width=width, height=height, save_figure=FALSE,
			significance_T1=significance_T1, significance_T2=significance_T2, loglog=loglog)
		p$p <- p$p + geom_label_repel(
			data = unique(subset(p$dt, y > transform_y(threshold)) %>% group_by(labels) %>% 
				filter(y == max(y))) %>% ungroup(),
			size = 2, aes(label=labels),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',
			min.segment.length=0)
		ggsave(filename=paste0("Figures/", file, '_categories.pdf'), p$p, width=width*scaling, height=height*scaling, units='mm')
		ggsave(filename=paste0("Figures/", file, '_categories.png'), p$p, width=width*scaling, height=height*scaling, units='mm')
	}
}

file_root <- c("meta_analysis")
for (i in 1:length(file_root)) {
	meta_list <- fread(
		paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
			file_root[i], "_figure_4.tsv.gz"))
	meta_list <- meta_list %>% filter(phenotype != "AlcCons")
	# Create the Burden, SKAT, and SKAT-O versions
	# Do the same thing, but split by case-control vs cts (way more power for cts).
	for (cc in c(TRUE, FALSE)) {
		cat(ifelse(cc, "case control\n", "continuous\n"))
		meta_list_tmp <- meta_list %>% filter(case_control == cc)
		for (phe in unique(meta_list_tmp$phenotype)) {
			cat(paste0(phe, "..."))
			meta_list_tmp_pheno <- meta_list_tmp %>% filter(phenotype == phe)
			p <- make_manhattan_plot(meta_list_tmp_pheno$chromosome_name,
				meta_list_tmp_pheno$start_position,
				meta_list_tmp_pheno$Pvalue,
				threshold=1000, significance_T = significance_T1,
				label=meta_list_tmp_pheno$external_gene_name, 
				colour_1 = "#6583E6",
				colour_2 = "#384980",
				loglog=TRUE)

			p$p <- p$p + geom_hline(yintercept=-log10(significance_T2), color='black', linetype='dashed')
			threshold <- ifelse(cc, 10, 50)

			if (nrow(p$dt %>% filter(y > transform_y(threshold))) < 20) {
				threshold <- 10
			}

			p$p <- p$p + geom_label_repel(
				data = unique(subset(p$dt, y > transform_y(threshold)) %>% group_by(labels) %>% 
					filter(y == max(y))) %>% ungroup(),
				size = 3, aes(label=labels),
				color='grey30', box.padding = 0.2, force = 0.3,
				label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',
				min.segment.length=0)
			width <- 230
			height <- 100
			scaling <- 1
			file <- paste0(file_root[i], "_", phe)
			ggsave(paste0("Figures/", file, '_manhattan.png'), p$p, width=width*scaling,
				height=height*scaling, units='mm')
		}
	}
}

# # This is for each phenotype
# for (i in 1:length(file_root)) {
# 	meta_list <- fread(
# 		paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
# 			file_root[i], "_figure_4.tsv.gz"))
# 	for (phe in unique(meta_list$phenotype)) {
# 		cat(phe, "\n")
# 		meta_list_tmp <- meta_list %>% filter(phenotype == phe)
# 		p <- make_manhattan_plot(meta_list_tmp$chromosome_name,
# 			meta_list_tmp$start_position,
# 			meta_list_tmp$Pvalue,
# 			threshold=1000, significance_T = 6.7e-7,
# 			label=meta_list_tmp$external_gene_name, 
# 			colour_1 = "#6583E6",
# 			colour_2 = "#384980")
# 		threshold <- ifelse(meta_list_tmp$case_control[1], 10, 10)
# 		p$p <- p$p + geom_label_repel(
# 			data = unique(subset(p$dt, y > threshold) %>% group_by(labels) %>% 
# 				filter(y == max(y))) %>% ungroup(),
# 			size = 2, aes(label=labels),
# 			color='grey30', box.padding = 0.2, force = 0.3,
# 			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50')
# 		width <- 230
# 		height <- 100
# 		scaling <- 1
# 		file <- paste0(phe, "_", file_root[i])
# 		ggsave(paste0("Figures/", file, '.jpg'), p$p, width=width*scaling,
# 			height=height*scaling, units='mm')
# 		width <- 150
# 		print(meta_list_tmp %>% filter(Pvalue < 6.7e-7))
# 	}
# }

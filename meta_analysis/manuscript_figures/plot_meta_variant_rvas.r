#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(biomaRt)
# devtools::install_github("mkanai/rgsutil")
library(rgsutil)
# BiocManager::install("GenomicRanges")
library(GenomicRanges)
source("../meta_analysis_utils.r")
source("../../QC/utils/pretty_plotting.r")

significance_T <- 8.0e-9
loglog <- TRUE

# gawk '$3 == "gene"' Homo_sapiens.GRCh38.114.gtf |
# gawk -F'\t' '{
#     match($9, /gene_id "([^"]+)"/, a);
#     match($9, /gene_name "([^"]+)"/, b);
#     print a[1] "\t" $1 "\t" $4 "\t" $5 "\t" b[1]
# }' > all_genes.tsv
genes_gtf_path <- "all_genes.tsv.gz"

loglog_trans <- function(y) {
	out <- numeric(length(y))
	log_part <- y <= 10
	out[log_part] <- y[log_part]
	out[!log_part] <- 10 + 2*(log2(y[!log_part] - 10 + 1))
	return(out)
}

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
		Musculoskeletal = color_musculoskeletal,
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

annotate_with_gene_names <- function(dt, significance_T=8e-9, return_binned=TRUE,
	biomaRt=FALSE, filter_to_genes=NULL)
{
	if (significance_T > 0) {
		significance_T <- log10(significance_T)
	}
	dt <- dt %>% filter(pval < significance_T)

	if (biomaRt) {
		# Connect to Ensembl
		mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror="www")

		# biomaRt wants chromosomes as "1", "2", ..., "X", "Y"
		genes <- getBM(
		  attributes = c("chromosome_name", "start_position", "end_position",
		                 "external_gene_name", "ensembl_gene_id"),
		  filters = c("chromosome_name", "start", "end"),
		  values = list(
		    chromosome_name = dt$contig,
		    start = dt$position,
		    end = dt$position
		  ),
		  mart = mart
		)

		genes$chr <- genes$chromosome_name
		genes$start <- genes$start_position
		genes$end <- genes$end_position

	} else {
		genes <- fread(genes_gtf_path) %>%
			dplyr::rename(
				ensembl_gene_id=V1,
				chr=V2,
				start=V3,
				end=V4,
				external_gene_name=V5)
	}

	genes <- unique(genes)

	if (!is.null(filter_to_genes)) {
		dt_filter <- filter_to_genes %>% dplyr::select(ensembl_gene_id)
		setkey(dt_filter, "ensembl_gene_id")
		setkey(genes, "ensembl_gene_id")
		genes <- merge(genes, dt_filter)
	}

	# Gene table (already queried from biomaRt)
	# Assume columns: chromosome_name, start_position, end_position,
	# ensembl_gene_id, external_gene_name
	
	# Make GRanges for variants
	variant_gr <- GRanges(
	  seqnames = dt$contig,
	  ranges = IRanges(
	  	start = dt$position,
	  	end = dt$position)
	)

	variant_gr <- unique(variant_gr)

	# Make GRanges for genes
	gene_gr <- GRanges(
	  seqnames = genes$chr,
	  ranges = IRanges(start = genes$start, end = genes$end),
	  ensembl_gene_id = genes$ensembl_gene_id,
	  external_gene_name = genes$external_gene_name
	)

	gene_gr <- unique(gene_gr)

	# Find overlaps
	hits <- findOverlaps(variant_gr, gene_gr)

	# Combine info
	variant_hits <- as.data.frame(variant_gr[queryHits(hits)])
	gene_hits <- as.data.frame(gene_gr[subjectHits(hits)])
	# Combine into one data frame
	result <- bind_cols(
	  variant_hits[, c("seqnames", "start")],
	  gene_hits[, c("ensembl_gene_id", "external_gene_name")]
	)

	# Rename for clarity
	colnames(result)[1:2] <- c("contig", "position")

	if (return_binned)
	{
		# Next, remove any genes that don't have names but that are 
		# overlap with other ensembl IDs that do have gene symbols

		result <- data.table(result %>% filter(external_gene_name!=""))
		setkeyv(result, c("contig", "position"))
		setkeyv(dt, c("contig", "position"))

		# Merge in the most significant association for each remaining set 
		# of variants that are all significant
		result <- data.table(merge(result, dt) %>% 
			group_by(phenotype_category, external_gene_name) %>% 
			slice_min(pval, with_ties = FALSE))

		result_binned <- result %>% mutate(
				contig = as.character(contig),
				bin = floor(position / 1e7) * 1e7) %>% 
			group_by(contig, bin) %>%
			slice_min(pval, n = 5, with_ties = FALSE) %>% ungroup()
		result_binned <- data.table(result_binned)

		return(result_binned)
	} else {
		result <- data.table(result)
		setkeyv(result, c("contig", "position"))
		setkeyv(dt, c("contig", "position"))
		result <- merge(result, dt, allow.cartesian=TRUE)
		setkeyv(result, c("ensembl_gene_id", "phenotype"))
		return(result)
	}
}

make_gene_manhattan_category_plot <- function(dt, buffer=100000000,
	chr_lengths=chr_lengths_38, significance_T=8.0e-9, ggplot_theme=theme_bw,
	save_figure=FALSE, file='file_out', scaling=1, width=230, height=100, 
	print_p=FALSE, loglog=TRUE, log_p_vals=TRUE)
{  	
    contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    start_end <- get_start_and_end(chr_lengths)
    dt_contigs <- data.frame(contig=contigs_, start=start_end$start,
    	end=start_end$end)

    # Now, we also want to shift by phenotype category
    if (!(all(c("Pvalue", "phenotype_category", "chromosome_name") %in%
    	names(dt)))) {
    	cat(paste0("Pvalue, phenotype_category, or chromosome_name not ", 
    		"present in the passed data.table\n"))
    }

    # Take the vector of phenotype categories
    # Order them alphabetically
    phenotype_categories <- sort(unique(dt$phenotype_category))
    phenotype_category_list <- list()
    i <- 1
   	for (phenotype_category in phenotype_categories) {
   		phenotype_category_list[[phenotype_category]] <- i
   		i <- i+1
   	}
    dt_categories <- data.frame(
    	phenotype_category = phenotype_categories,
    	start = ((dt_contigs$end[nrow(dt_contigs)]) * 
    		(seq(1,length(phenotype_categories))-1) + 1),
    	end = ((dt_contigs$end[nrow(dt_contigs)]) * 
    		(seq(1,length(phenotype_categories))))) %>%
    	mutate(middle = floor(start + (end-start)/2), length = (end-start)) %>% 
    	mutate(shifted_position = middle + 
    		(unlist(phenotype_category_list[phenotype_category])-1) * buffer)

    dt_plot <- dt %>% transmute(
    	contig=gsub("chr", "", chromosome_name),
    	position=as.integer(position),
    	pval=as.numeric(Pvalue),
    	# labels=external_gene_name,
    	phenotype_category=phenotype_category) %>% mutate(
    	x = dt_categories[
    		unlist(phenotype_category_list[phenotype_category]), 'start'] + 
    		dt_contigs[gsub('X', 23, contig), 'start'] + 
    		position + 
    		(unlist(phenotype_category_list[phenotype_category])-1) * buffer)

    if (!log_p_vals) {
	  	dt_plot <- dt_plot %>% mutate(
	  		y = ifelse(pval < 1e-300, 300, -log10(pval)))
	} else {
		dt_plot <- dt_plot %>% mutate(y = -pval)
	}

  	if (loglog) {
		transform_y <- loglog_trans
	} else {
		transform_y <- function(y) { return(y) }
	}

	# Apply transformation
	dt_plot <- dt_plot %>% mutate(y_trans = transform_y(y))
	breaks <- c(2, 4, 6, 8, 10, 20, 50, 100, 200, 400, 800)
	breaks_trans <- transform_y(breaks)

    p <- ggplot(dt_plot, aes(x=x, y=y_trans, col=phenotype_category)) + 
    	geom_point_rast(size=0.5)

    if (loglog) {
	    p <- p + scale_y_continuous(
	    	breaks = breaks_trans,
	    	labels = breaks
	    )
    } else {
    	p <- p + scale_y_continuous(breaks=scales::pretty_breaks(n=10))
    }
    p <- p + geom_hline(yintercept=-log10(significance_T), color='#E15759',
    	linetype='dashed') +
        scale_x_continuous(
        	breaks=dt_categories$shifted_position,
        	labels=dt_categories$phenotype_category) +
        scale_color_manual(values = get_category_colors()) +
        labs(x='', y=TeX("$-\\log_{10}(P)$")) + ggplot_theme() + 
        theme(
        	axis.text.x = element_text(angle = 45, hjust = 1),
        	legend.position="none")

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling,
        	height=height*scaling, units='mm')
    }

    if (print_p) { print(p) }

    return(list(p=p, dt=dt_plot))
}

# BMRC
meta_list <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/manhattan_plots_null_results_downsampled.tsv.gz") %>% 
	dplyr::rename(chr=CHR, pos=POS, ref=REF, alt=ALT) %>% 
	mutate(chr=gsub("chr", "", chr)) %>% 
	mutate(phenotype_category = unlist(phenotype_broad_categories[phenotype])) %>% 
	mutate(case_control = ifelse(type == "binary", TRUE, FALSE))

# Do the same thing, but split by case-control vs cts (way more power for cts).
for (anc in unique(meta_list$ancestry)) {
	for (cc in c(TRUE, FALSE)) {
		cat(ifelse(cc, "case control\n", "continuous\n"))
		meta_list_tmp <- meta_list %>% filter(case_control == cc, ancestry == anc)
		p <- make_manhattan_plot(meta_list_tmp$chr,
			meta_list_tmp$pos,
			-meta_list_tmp$`P-value`,
			threshold=1000,
			significance_T = significance_T,
			label=NULL, 
			colour_1 = "#6583E6",
			colour_2 = "#384980",
			loglog=!cc,
			log_p_vals=TRUE)

		dt_label <- annotate_with_gene_names(meta_list_tmp %>%
			dplyr::rename(contig = chr, position=pos, pval=`P-value`)) %>% 
			mutate(y=-pval)
		if (!cc) {
			dt_label$y <- loglog_trans(dt_label$y)
		}
		dt_label <- data.table(dt_label %>% group_by(ID) %>% 
			slice_min(pval, with_ties = FALSE) %>% ungroup())

		p$dt <- data.table(p$dt)
		setkeyv(p$dt, c("contig", "position", "y"))
		setkeyv(dt_label, c("contig", "position", "y"))
		dt_label <- merge(p$dt, dt_label)

		p$p <- p$p + geom_label_repel(
			data = dt_label,
			size = 3, aes(label=external_gene_name),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',
			min.segment.length=0)
		width <- 230
		height <- 100
		scaling <- 1
		file <- paste0("meta_analysis", "_",
			ifelse(cc, "variant_case_control", "variant_cts"))
		ggsave(paste0("Figures/", file,  '_', anc, '.pdf'), p$p, width=width*scaling,
			height=height*scaling, units='mm')
		print("done creating!")

		width <- 150
		p <- make_gene_manhattan_category_plot(
			meta_list_tmp %>% dplyr::rename(
				position = pos, Pvalue=`P-value`, chromosome_name=chr),
			buffer=1000000000,
			scaling=scaling, width=width, height=height, save_figure=FALSE,
			significance_T=significance_T, loglog=!cc)

		# Could consider including genes
		dt_label <- annotate_with_gene_names(p$dt) %>% 
			mutate(y=-pval) %>% dplyr::select(-pval)
		p$dt <- data.table(p$dt)
		setkeyv(p$dt, intersect(names(p$dt), names(dt_label)))
		setkeyv(dt_label, intersect(names(p$dt), names(dt_label)))
		dt_label <- merge(p$dt, dt_label)

		p$p <- p$p + geom_label_repel(
			data = dt_label,
			size = 3, aes(label=external_gene_name),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50',
			min.segment.length=0)

		ggsave(
			filename=paste0("Figures/", file, '_', anc, '_categories.pdf'), p$p,
			width=width*scaling, height=height*scaling, units='mm')
	}
}

# # Plot the results against each other (where possible)
# # Note that we need to also include the remaining cohorts in the
# # rare-variant meta-analysis.

# # BMRC
# # dt_gene <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_all_meta_subsets_060625.tsv.gz")
# # dt_gene <- dt_gene %>% mutate(phenotype = gsub("_ALL", "", phenotype))
# dt_gene <- fread(
# 	"significant_assocs_from_all_meta_subsets_060625.tsv.gz") %>%
# 	mutate(phenotype = gsub("_ALL", "", phenotype))
# dt_gene <- fread(
# 	"all_assocs_from_meta_060625.tsv.gz") %>%
# 	mutate(phenotype = gsub("_ALL$", "", phenotype))
# # BMRC
# # meta_list <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/variant/manhattan_plots.tsv.gz")
# # meta_list <- meta_list %>% 
# # 	dplyr::rename(chr=CHR, pos=POS, ref=REF, alt=ALT) %>% 
# # 	mutate(chr=gsub("chr", "", chr))
# dt_variant <- fread("manhattan_plots.tsv.gz") %>% 
# 	dplyr::rename(chr=CHR, pos=POS, ref=REF, alt=ALT) %>% 
# 	dplyr::mutate(
# 		chr=gsub("chr", "", chr),
# 		phenotype_category = unlist(phenotype_broad_categories[phenotype]),
# 		case_control = ifelse(type == "binary", TRUE, FALSE)) %>%
# 	dplyr::rename(contig = chr, position=pos, pval=`P-value`)

# dt_variant <- annotate_with_gene_names(dt_variant, significance_T=1,
# 	return_binned=FALSE,
# 	filter_to_genes=unique(
# 		dt_gene %>% filter(ancestry == "all") %>%
# 		dplyr::rename(ensembl_gene_id=Region) %>% dplyr::select(ensembl_gene_id))
# 	)

# # Now we can plot them against each other
# dt_gene <- dt_gene %>% filter(ancestry == "all") %>% 
# 	mutate(pval_gene=log10(Pvalue)) %>% dplyr::select(-Pvalue)
# dt_variant <- dt_variant %>% mutate(
# 	pval_variant=pval,
# 	Region=ensembl_gene_id) %>% dplyr::select(-c("pval", "ensembl_gene_id"))

# dt_gene <- dt_gene %>% group_by(Region, phenotype) %>% slice_min(pval_gene)
# dt_variant <- dt_variant %>% group_by(Region, phenotype) %>% slice_min(pval_variant)
# dt_gene <- data.table(dt_gene)
# dt_variant <- data.table(dt_variant)

# setkeyv(dt_gene, c("Region", "phenotype"))
# setkeyv(dt_variant, c("Region", "phenotype"))

# dt_plot <- merge(dt_gene, dt_variant)
# ggplot_theme <- theme_bw

# for (cc in c(TRUE, FALSE)) {
# 	cat(ifelse(cc, "case control\n", "continuous\n"))

# 	p <- ggplot(dt_plot %>% filter(case_control == cc),
# 		aes(x=-pval_variant, y=-pval_gene, col=phenotype_category)) + 
# 		geom_point_rast(size=0.5) + 
# 		scale_color_manual(values = get_category_colors()) +
# 		labs(x=TeX("$-\\log_{10}(P_{variant})$"),
# 			y=TeX("$-\\log_{10}(P_{gene})$")) + ggplot_theme() +
# 		geom_abline(
# 			slope = 1,
# 			intercept = 0,
# 			linetype = "dashed",
# 			color = "gray") +
# 		scale_x_continuous(limits = c(0, NA)) +
#  		scale_y_continuous(limits = c(0, NA)) + 
#  		theme(legend.title=element_blank())
# 	# theme(legend.position="none")
# 	# plot(-dt_plot$pval_variant, -dt_plot$pval_gene)
# 	# abline(0,1, col='red')
# 	print(p)

# 	p <- ggplot(dt_plot %>% filter(case_control == cc),
# 		aes(x=-pval_variant, y=-pval_gene, col=phenotype_category)) + 
# 		geom_point_rast(size=0.5) + 
# 		scale_color_manual(values = get_category_colors()) +
# 		labs(x=TeX("$-\\log_{10}(P_{variant})$"),
# 			y=TeX("$-\\log_{10}(P_{gene})$")) + ggplot_theme() +
# 		geom_abline(
# 			slope = 1,
# 			intercept = 0,
# 			linetype = "dashed",
# 			color = "gray") +
# 		scale_x_continuous(limits = c(0, 20)) +
#  		scale_y_continuous(limits = c(0, 20)) + 
#  		theme(legend.title=element_blank())
#  	print(p)
# }

# # CHECK THIS FOR MAF < 0.1%

# # THIS SUBSEQUENT PART IS FOR COUNTING EFFECTS
# # 	dt_sig <- data.table(dt_sig)
# # 	setkeyv(dt_sig, c("chr", "pos", "ref", "alt"))

# # 	result <- data.table(result)
# # 	setkeyv(result, c("chr", "pos", "ref", "alt"))

# # 	dt_final <- merge(dt_sig, result)
# # 	dt_final <- dt_final %>% dplyr::rename(Region = "ensembl_gene_id")
# # 	dt_final <- data.table(dt_final)
# # 	setkeyv(dt_final, c("Region", "phenotype"))

# # # Next, merge back in the traits that are associated
# # dt_gene <- fread(
# # 	"../significant_assocs_from_all_meta_subsets_060625.tsv.gz") %>%
# # 	mutate(phenotype = gsub("_ALL", "", phenotype))

# # # We now use the gene, phenotypes define the match
# # dt_gene <- dt_gene %>% filter(ancestry == "all")
# # dt_gene <- unique(data.table(dt_gene %>% dplyr::select(Region, phenotype)))
# # setkeyv(dt_gene, c("Region", "phenotype"))
# # dt_overlap <- merge(dt_gene, dt_final)
# # # And then determine which (gene, phenotype) pairs have...
# # # A single variant association driving the signal
# # # More than one
# # # None

# # # Finally, merge in the gene-based associations that are significiant 
# # # in the meta-analysis, and see what the overlap is.

# # This is to generate all comparisons.
# # manhattan_plots_to_merge.tsv.gz is dt before the merge at the end of 
# # annotate_with_gene_names
# # variant_gene_map.tsv.gz is result before the merge at the end of 
# # annotate_with_gene_names. We need to split by chromosome below to avoid 
# # running out of memory. Really this shouldn't be done with R, but easier than
# # rewriting.

# # DEV: the following is not quite right, because we have not restricted to the collection
# # of variants that are present in each of the group based tests.

# dt <- fread("manhattan_plots_to_merge.tsv.gz")
# setkeyv(dt, c("contig", "position"))
# result <- fread("variant_gene_map.tsv.gz")
# setkeyv(result, c("contig", "position"))
# dt_gene <- fread(
# 	"all_assocs_from_meta_060625.tsv.gz") %>%
# 	mutate(phenotype = gsub("_ALL$", "", phenotype))
# 	dt_gene <- dt_gene %>% filter(ancestry == "all") %>% 
# 	mutate(pval_gene=log10(Pvalue)) %>% dplyr::select(-Pvalue)
# dt_gene <- dt_gene %>% group_by(Region, phenotype) %>% slice_min(pval_gene)
# dt_gene <- data.table(dt_gene)
# setkeyv(dt_gene, c("Region", "phenotype"))

# dt_plot_list <- list()
# for (chr in as.character(c(seq(1,22), "X"))) {
# 	cat(paste(chr, "\n"))
# 	dt_chr <- dt[contig == chr]
# 	dt_chr <- merge(dt_chr, result)
# 	setkeyv(dt_chr, c("ensembl_gene_id", "phenotype"))
# 	dt_chr <- dt_chr %>% mutate(
# 		pval_variant=pval,
# 		Region=ensembl_gene_id) %>% dplyr::select(-c("pval", "ensembl_gene_id"))
# 	dt_chr <- dt_chr %>% group_by(Region, phenotype) %>% slice_min(pval_variant)
# 	dt_chr <- data.table(dt_chr)
# 	setkeyv(dt_chr, c("Region", "phenotype"))
# 	print(nrow(dt_chr))
# 	dt_plot_list[[chr]] <- merge(dt_gene, dt_chr)
# 	print(nrow(dt_plot_list[[chr]]))
# }

# dt_plot <- rbindlist(dt_plot_list)

# ggplot_theme <- theme_bw
# for (cc in c(TRUE, FALSE)) {
# 	cat(ifelse(cc, "case control\n", "continuous\n"))

# 	p <- ggplot(dt_plot %>% filter(case_control == cc),
# 		aes(x=-pval_variant, y=-pval_gene, col=phenotype_category)) + 
# 		geom_point_rast(size=0.5) + 
# 		scale_color_manual(values = get_category_colors()) +
# 		labs(x=TeX("$-\\log_{10}(P_{variant})$"),
# 			y=TeX("$-\\log_{10}(P_{gene})$")) + ggplot_theme() +
# 		geom_abline(
# 			slope = 1,
# 			intercept = 0,
# 			linetype = "dashed",
# 			color = "gray") +
# 		scale_x_continuous(limits = c(0, NA)) +
#  		scale_y_continuous(limits = c(0, NA)) + 
#  		theme(legend.title=element_blank())
# 	# theme(legend.position="none")
# 	# plot(-dt_plot$pval_variant, -dt_plot$pval_gene)
# 	# abline(0,1, col='red')
# 	print(p)

# 	p <- ggplot(dt_plot %>% filter(case_control == cc),
# 		aes(x=-pval_variant, y=-pval_gene, col=phenotype_category)) + 
# 		geom_point_rast(size=0.5) + 
# 		scale_color_manual(values = get_category_colors()) +
# 		labs(x=TeX("$-\\log_{10}(P_{variant})$"),
# 			y=TeX("$-\\log_{10}(P_{gene})$")) + ggplot_theme() +
# 		geom_abline(
# 			slope = 1,
# 			intercept = 0,
# 			linetype = "dashed",
# 			color = "gray") +
# 		scale_x_continuous(limits = c(0, 20)) +
#  		scale_y_continuous(limits = c(0, 20)) + 
#  		theme(legend.title=element_blank())
#  	print(p)
# }


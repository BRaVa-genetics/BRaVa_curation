#!/bin/Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(future.apply)
library(progressr)
# devtools::install_github("mkanai/rgsutil")
# library(rgsutil)
source("../meta_analysis_utils.r")
source("../../QC/utils/pretty_plotting.r")

significance_T1 <- 0.05/(20000*3*3*2)
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

make_gene_manhattan_category_plot <- function(dt, buffer=100000000,
	chr_lengths=chr_lengths_38, significance_T1=6.7e-7, significance_T2=NULL, ggplot_theme=theme_bw,
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
    p <- p + geom_hline(yintercept=-log10(significance_T1), color='#E15759', linetype='dashed') +
    	# geom_hline(yintercept=-log10(significance_T2), color='black', linetype='dashed') +
        scale_x_continuous(breaks=dt_categories$shifted_position, labels=dt_categories$phenotype_category) +
        scale_color_manual(values = get_category_colors()) +
        labs(x='', y = expression(-log[10](italic(P)))) + ggplot_theme() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")

    if (save_figure) {
        ggsave(paste0(file, '.png'), p, width=width*scaling, height=height*scaling, units='mm')
        ggsave(paste0(file, '.pdf'), p, width=width*scaling, height=height*scaling, units='mm')
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

plan(multisession, workers = 16)   # adjust as needed
handlers("progress")

project_dir <- getwd()

file_root <- c("meta_analysis", "AFR", "AMR", "EAS", "EUR", "SAS", "non_EUR")

job_dt <- CJ(
  file_root = file_root,
)

with_progress({
  p <- progressor(along = seq_len(nrow(job_dt)))

  future_lapply(seq_len(nrow(job_dt)), function(i) {

    setwd(project_dir)

    root <- job_dt$file_root[i]

    p(sprintf("%s", root))

    meta_list <- fread(
      paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
             root, "_figure_4.tsv.gz")
    )
    setDT(meta_list)

    meta_list[, phenotype_category := unlist(phenotype_broad_categories[phenotype])]
    meta_list[, external_gene_name := fifelse(
      is.na(external_gene_name),
      ensembl_gene_id,
      external_gene_name
    )]

    for (cc in c(TRUE, FALSE)) {

      cat(ifelse(cc, "case control\n", "continuous\n"))

      meta_list_tmp <- meta_list[case_control == cc]
      if (nrow(meta_list_tmp) == 0L) next

      p1 <- make_manhattan_plot(
        meta_list_tmp$chromosome_name,
        meta_list_tmp$start_position,
        meta_list_tmp$Pvalue,
        threshold = 1000,
        significance_T = significance_T1,
        label = meta_list_tmp$external_gene_name,
        colour_1 = "#6583E6",
        colour_2 = "#384980",
        loglog = TRUE
      )

      dt_plot1 <- as.data.table(p1$dt)

      threshold <- ifelse(cc, 10, 50)

      if (nrow(dt_plot1[y > transform_y(threshold)]) < 20) {
        threshold <- 10
      }

      gene_labels_to_plot <- unique(
        dt_plot1[y > transform_y(threshold)]
      )[
        , .SD[which.max(y)], by = labels
      ]

      p1$p <- p1$p + geom_label_repel(
        data = gene_labels_to_plot,
        aes(label = labels),
        fontface = "italic",
        size = 3,
        color = "grey30",
        box.padding = 0.2,
        force = 1,
        label.padding = 0.1,
        point.padding = 0.1,
        segment.color = "grey50",
        min.segment.length = 0,
        segment.size = 0.5,
        segment.alpha = 0.8
      )

      width <- 230
      height <- 100
      scaling <- 1

      file <- paste0(root, "_", ifelse(cc, "unique_case_control", "unique_cts"))

      cat(paste0(file, "\n"))

      ggsave(
        paste0("Figures/", file, ".pdf"),
        p1$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )
      ggsave(
        paste0("Figures/", file, ".png"),
        p1$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )

      width <- 150

      p2 <- make_gene_manhattan_category_plot(
        meta_list_tmp,
        buffer = 1000000000,
        scaling = scaling,
        width = width,
        height = height,
        save_figure = FALSE,
        significance_T1 = significance_T1,
        significance_T2 = significance_T2,
        loglog = loglog
      )

      dt_plot2 <- as.data.table(p2$dt)

      gene_labels_to_plot2 <- unique(
        dt_plot2[y > threshold]
      )[
        , .SD[which.max(y)], by = .(labels, phenotype_category)
      ]

      p2$p <- p2$p + geom_label_repel(
        data = gene_labels_to_plot2,
        aes(label = labels),
        size = 2,
        fontface = "italic",
        color = "grey30",
        box.padding = 0.2,
        force = 1,
        label.padding = 0.1,
        point.padding = 0.1,
        segment.color = "grey50",
        min.segment.length = 0
      )

      ggsave(
        filename = paste0("Figures/", file, "_categories.png"),
        p2$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )
      ggsave(
        filename = paste0("Figures/", file, "_categories.pdf"),
        p2$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )
    }

    NULL
  }, future.seed = TRUE)
})

plan(multisession, workers = 16)   # adjust as needed
handlers("progress")

file_root <- c("meta_analysis")

for (root in file_root) {

  meta_list <- fread(
    paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/",
           root, "_figure_4.tsv.gz")
  )
  setDT(meta_list)

  meta_list[, external_gene_name := fifelse(
    is.na(external_gene_name),
    ensembl_gene_id,
    external_gene_name
  )]

  phenotypes <- unique(meta_list$phenotype)

  with_progress({
    p <- progressor(along = phenotypes)

    future_lapply(phenotypes, function(phe) {

      setwd(project_dir)

      p(phe)

      meta_list_tmp <- meta_list[phenotype == phe]
      if (nrow(meta_list_tmp) == 0L) return(NULL)

      p <- make_manhattan_plot(
        meta_list_tmp$chromosome_name,
        meta_list_tmp$start_position,
        meta_list_tmp$Pvalue,
        threshold = 1000,
        significance_T = significance_T1,
        label = meta_list_tmp$external_gene_name,
        colour_1 = "#6583E6",
        colour_2 = "#384980",
        loglog = TRUE
      )

      p$dt <- as.data.table(p$dt)
      p$dt[, y := ifelse(y > 300, 300, y)]

      threshold <- -log10(significance_T1)

      if (any(p$dt$y > transform_y(threshold))) {
        gene_labels_to_plot <- p$dt[
          y > transform_y(threshold)
        ][
          , .SD[which.max(y)], by = labels
        ]

        cat(paste0("number of significant genes: ", nrow(gene_labels_to_plot), "\n"))

        if (nrow(gene_labels_to_plot) > 50) {
          gene_labels_to_plot <- gene_labels_to_plot[order(-y)][1:50]
        }

        p$p <- p$p + geom_label_repel(
          data = gene_labels_to_plot,
          aes(label = labels),
          fontface = "italic",
          size = 3,
          color = "grey30",
          label.padding = 0.1,
          box.padding = 0.2,
          point.padding = 0.1,
          force = 1,
          max.overlaps = Inf,
          segment.color = "grey50",
          segment.size = 0.5,
          segment.alpha = 0.8,
          min.segment.length = 0
        )
      }

      width <- 150
      height <- 80
      scaling <- 1
      file <- paste0(phe, "_", root)

      ggsave(
        paste0("Figures/", file, ".png"),
        p$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )
      ggsave(
        paste0("Figures/", file, ".pdf"),
        p$p,
        width = width * scaling,
        height = height * scaling,
        units = "mm"
      )

      NULL
    }, future.seed = TRUE)
  })
}

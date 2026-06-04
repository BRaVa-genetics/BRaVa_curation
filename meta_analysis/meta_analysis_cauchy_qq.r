#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)

source("../QC/utils/pretty_plotting.r")
gene_name_mapping_file <- "../data/gene_mapping.txt.gz"
T <- 6 # -log10(P) threshold for inclusion on the plots

main <- function(args)
{
    ribbon_p <- 0.95
    cex_labels <- 2

    if (!is.null(args$meta_analysis_results_folder)) {
        files <- dir(args$meta_analysis_results_folder, full.names=TRUE)
        files_to_include <- grep("meta", dir(args$meta_analysis_results_folder))
        files <- files[files_to_include]
    } else if (!is.null(args$meta_analysis_results_file)) {
        files <- args$meta_analysis_results_file
    } else {
        stop("Neither a folder or single file of results have been provided to plot")
    }

    if (args$include_gene_names) {
        if (file.exists(gene_name_mapping_file)) {
            dt_genes <- fread(gene_name_mapping_file)
            if (!(all(c("Gene stable ID", "Gene name") %in% names(dt_genes)))) {
                warning(paste0("gene name mapping file does not contain 'Gene stable ID' ",
                    "(ensembl ID *without* version) and 'Gene name' (gene symbol for plotting), ",
                    "reverting to not plotting gene-names"))
                 args$include_gene_names <- FALSE
            } else {
                names(dt_genes)[which(names(dt_genes) == "Gene stable ID")] <- "Region"
                names(dt_genes)[which(names(dt_genes) == "Gene name")] <- "labels"
                setkey(dt_genes, "Region")
            }
        } else {
            warning("gene name mapping file is not present, reverting to not plotting gene-names")
            args$include_gene_names <- FALSE
        }
    }

    width_plot <- ifelse(is.null(args$type), 8, 4)
    width_plot <- ifelse(args$burden_only_plot, 6, width_plot)

    pdf(file=args$out, width=width_plot, height=4)
    for (file in files) {
        phe <- gsub(".*/(.*)_meta.*", "\\1", file)
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", phe))))
        cat(paste0(phe_plot, "...\n"))
        dt_meta <- fread(file, key="Region")
        dt_meta[, Pvalue := ifelse(Pvalue == 0, 320, -log10(Pvalue))]
        if(args$include_gene_names) {
            dt_meta <- merge(dt_meta, dt_genes, all.x=TRUE)
        }
        
        if (!is.null(args$type)) {
            dt_meta <- dt_meta %>% filter(type == args$type)
        }

        dt_meta_to_plot <- data.table(
            dt_meta %>% 
            arrange(class, type, Group, desc(Pvalue)) %>%
            select(-c(class, type, Group)),
            dt_meta %>% 
            group_by(class, type, Group) %>% 
            arrange(desc(Pvalue)) %>% 
            summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
        	)
        )

        for (g in unique(dt_meta_to_plot$Group))
        {
            p <- create_pretty_qq_plot(
                plot_title=phe_plot,
                plot_subtitle=paste0("Cauchy: ", g),
                cex_labels=cex_labels,
                dt_meta_to_plot %>% filter(Group == g),
                aes(x=Pvalue_expected, y=Pvalue, color=class),
                save_figure=FALSE,
                x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                y_label=TeX("$-\\log_{10}(P_{observed})$"),
                key_cols=c("class", "Pvalue"),
                aes_ribbon = aes(ymin=clower, ymax=cupper),
                print_p=FALSE)

            if(args$include_gene_names) {
                p <- p + geom_label_repel(data=dt_meta_to_plot %>% 
                    filter(Group == g, class == "Burden", Pvalue > T),
                    aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                    color = 'grey30', segment.color = 'grey50',
                    size=cex_labels, segment.size=0.1, show.legend = FALSE)
            } 
            p <- p + facet_wrap(~type)
            print(p)
        }
    }
    dev.off()
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--meta_analysis_results_folder",
    default=NULL, required=FALSE,
    help="Folder containing the meta-analysis results files")
parser$add_argument("--meta_analysis_results_file",
    default=NULL, required=FALSE,
    help="Meta-analysis results file")
parser$add_argument("--out",
    default="meta_analysis_qq_100.pdf",
    required=FALSE, help="Output filepath")
parser$add_argument("--type", default=NULL, required=FALSE,
    help="Which meta-analysis results to plot {'Stouffer', 'Weighted Fisher', 'Inverse variance weighted'}")
parser$add_argument("--include_gene_names", default=FALSE, action='store_true',
    help="Do you want to highlight the most significant genes with their gene-names?")
parser$add_argument("--burden_only_plot", default=FALSE, action='store_true',
    help="Do you want to create plots containing just the burden p-values together with colour-coded effect direction?")
args <- parser$parse_args()

main(args)

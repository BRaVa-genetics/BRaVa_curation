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

source("meta_analysis_utils.r")

main <- function(args)
{
    print(args)
    ribbon_p <- 0.95

    if (!is.null(args$analysis_results_folder)) {
        files <- dir(args$analysis_results_folder,
            pattern=args$analysis_results_path_regexp,
            full.names=TRUE)
    } else if (!is.null(args$analysis_results_file)) {
        files <- args$analysis_results_file
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

    pdf(file=args$out, width=6, height=4)
    for (file in files) {
        cat(paste0(file, "\n"))
        file_info <- extract_file_info(gsub(".*/(.*)", "\\1", file))
        phe_plot <- gsub("_", " ", gsub("_$", "", str_trim(gsub("[[:space:]_]+", "\\_", file_info$phenotype))))
        phe_plot <- paste(file_info$dataset, file_info$ancestry, phe_plot, sep=", ")
        if (file_info$binary) {
            phe_plot <- paste0(phe_plot,
                "\nn cases = ", prettyNum(file_info$n_cases, big.mark=",", trim=TRUE),
                ", n controls = ", prettyNum(file_info$n_controls, big.mark=",", trim=TRUE)
            )
        } else {
            phe_plot <- paste0(phe_plot, "\nn = ", prettyNum(file_info$n, big.mark=",", trim=TRUE))
        }
        cat(paste0(phe_plot, "...\n"))
        dt <- fread(file)
        if (args$burden_only_plot) {
            dt <- dt %>% select(-c("Pvalue", "Pvalue_SKAT")) %>% 
                rename(Pvalue = Pvalue_Burden) %>%
                mutate(Pvalue = ifelse(Pvalue == 0, 320, -log10(Pvalue)), class = "Burden") %>% 
                select(Region, Group, max_MAF, Pvalue, BETA_Burden, SE_Burden, class)
        } else {
            dt <- melt(dt, id.vars = c("Region", "Group", "max_MAF"),
                    measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
                    value.name = "Pvalue", variable.name = "class") %>% 
            mutate(Pvalue = ifelse(Pvalue == 0, 320, -log10(Pvalue)))
        }
        dt <- data.table(dt)
        setkey(dt, "Region")

        if(args$include_gene_names) {
            dt <- merge(dt, dt_genes, all.x=TRUE)
        }
        
        dt[, class := as.factor(class)]
        dt$max_MAF[which(is.na(dt$max_MAF))] <- "Cauchy"
        levels(dt$class) <- list(
            `SKAT-O` = "Pvalue",
            Burden = "Pvalue_Burden",
            SKAT = "Pvalue_SKAT")

        dt_to_plot <- data.table(
            dt %>% 
            arrange(Group, max_MAF, class, desc(Pvalue)) %>%
            select(-c(Group, max_MAF, class)),
            dt %>% 
            group_by(Group, max_MAF, class) %>% 
            arrange(desc(Pvalue)) %>% 
            summarize(
                Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
                clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
                cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
        	)
        )

        if (args$burden_only_plot) {
            dt_to_plot <- dt_to_plot %>% mutate(OR = exp(BETA_Burden))
            dt_to_plot$color <- cut(dt_to_plot$OR,
                breaks = c(-Inf, 0.95, 1, 1.05, Inf),
                labels = c("< 0.95", "[0.95, 1)", "[1, 1.05]", "> 1.05"))
            dt_to_plot$color <- factor(dt_to_plot$color,
                levels = c("< 0.95", "[0.95, 1)", "[1, 1.05]", "> 1.05"))
            dummy_data <- data.frame(Pvalue_expected = NA, Pvalue = NA,
                color = factor(c("< 0.95", "[0.95, 1)", "[1, 1.05]", "> 1.05"),
                levels = c("< 0.95", "[0.95, 1)", "[1, 1.05]", "> 1.05"))
            )
        }

        max_MAF_groups <- setdiff(unique(dt_to_plot$max_MAF), "Cauchy")
        groups <- setdiff(unique(dt_to_plot$Group), "Cauchy")
        for (m in max_MAF_groups) {
            for (g in groups) {
                # max_MAF_plot <- ifelse(
                #     grepl("e", as.character(m)),
                #     gsub("([0-9\\.]+)e(-)*0*([1-9][0-9]*)", "$\\1\\\\times 10^{\\2\\3}$", as.character(m)),
                #     paste0("$", as.character(m), "$"))
                max_MAF_plot <- as.character(m)
                variant_class_plot <- gsub("_", " ", gsub("[\\|;]", ",\n", g))
                cex_labels <- 2

                if (!args$burden_only_plot)
                {
                    p <- create_pretty_qq_plot(
                        plot_title=phe_plot,
                        plot_subtitle=paste0(variant_class_plot, "; max MAF = ", max_MAF_plot),
                        cex_labels=cex_labels,
                        dt_to_plot %>% filter(Group==g, max_MAF==m),
                        aes(x=Pvalue_expected, y=Pvalue, color=class),
                        save_figure=FALSE,
                        x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                        y_label=TeX("$-\\log_{10}(P_{observed})$"),
                        key_cols=c("class", "Pvalue"),
                        aes_ribbon = aes(ymin=clower, ymax=cupper),
                        print_p=FALSE
                    )

                    if(args$include_gene_names) {
                        p <- p + geom_label_repel(data=dt_to_plot %>% 
                            filter(Group==g, max_MAF==m, class == "Burden", Pvalue > T),
                            aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                            color = 'grey30', segment.color = 'grey50',
                            size=cex_labels, segment.size=0.1, show.legend = FALSE)
                    }
                    print(p)

                } else {
                    p <- create_pretty_qq_plot(
                        plot_title=phe_plot,
                        plot_subtitle=paste0(variant_class_plot, "; max MAF = ", max_MAF_plot),
                        cex_labels=cex_labels,
                        rbind(dt_to_plot %>% filter(Group==g, max_MAF==m, class=="Burden"), dummy_data, fill=TRUE),
                        aes(x=Pvalue_expected, y=Pvalue, color=color),#, size=size),
                        save_figure=FALSE,
                        x_label=TeX("$-\\log_{10}(P_{expected})$"), 
                        y_label=TeX("$-\\log_{10}(P_{observed})$"),
                        key_cols=c("class", "Pvalue"),
                        aes_ribbon = aes(ymin=clower, ymax=cupper),
                        print_p=FALSE
                    )

                    p <- p + scale_color_manual(
                        values = c(
                            "< 0.95" = "blue3",
                            "[0.95, 1)" = "cornflowerblue",
                            "[1, 1.05]" = "indianred3",
                            "> 1.05" = "red"),
                        labels = c("< 0.95" = "< 0.95",
                            "[0.95, 1)" = "[0.95, 1)",
                            "[1, 1.05]" = "[1, 1.05]",
                            "> 1.05" = "> 1.05"),
                        name = "Odds ratio",  aesthetics = c("colour", "fill")
                    ) + guides(colour = guide_legend(override.aes = list(size=5)))
                   
                    if(args$include_gene_names) {
                        p <- p + geom_label_repel(data=dt_to_plot %>% 
                            filter(Group==g, max_MAF==m, class == "Burden", Pvalue > T),
                            aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                            color = 'grey30', segment.color = 'grey50',
                            size=cex_labels, segment.size=0.1, show.legend = FALSE)
                    }
                    print(p)
                }
            }
        }

        p <- create_pretty_qq_plot(
            plot_title=phe_plot,
            plot_subtitle="Cauchy",
            cex_labels=2,
            dt_to_plot %>% filter(Group == "Cauchy"),
            aes(x=Pvalue_expected, y=Pvalue, color=class),
            save_figure=FALSE,
            x_label=TeX("$-\\log_{10}(P_{expected})$"), 
            y_label=TeX("$-\\log_{10}(P_{observed})$"),
            key_cols=c("class", "Pvalue"),
            aes_ribbon = aes(ymin=clower, ymax=cupper),
            print_p=FALSE)

        if(args$include_gene_names) {
            p <- p + geom_label_repel(data=dt_to_plot %>% 
                filter(Group=="Cauchy", class == "Burden", Pvalue > T),
                aes(label=labels), box.padding = 0.5, label.padding=0.1, point.padding = 0.2,
                color = 'grey30', segment.color = 'grey50',
                size=cex_labels, segment.size=0.1, show.legend = FALSE)
        }
        print(p)
    }
    dev.off()
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--analysis_results_folder",
    default=NULL, required=FALSE,
    help="Folder containing the analysis results files")
parser$add_argument("--analysis_results_file",
    default=NULL, required=FALSE,
    help="analysis results file")
parser$add_argument("--analysis_results_path_regexp", default="*cleaned*", required=FALSE,
    help="Regular expression of analysis results files to loop over")
parser$add_argument("--out", default="biobank_qq.pdf", required=FALSE,
    help="Output filepath")
parser$add_argument("--include_gene_names", default=FALSE, action='store_true',
    help="Do you want to highlight the most significant genes with their gene-names?")
parser$add_argument("--burden_only_plot", default=FALSE, action='store_true',
    help="Do you want to create plots containing just the burden p-values together with colour-coded effect direction?")
args <- parser$parse_args()

main(args)

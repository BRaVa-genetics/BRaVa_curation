#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(latex2exp)
library(argparse)
library(ggthemes)
library(scales)

source("../QC/utils/pretty_plotting.r")
source("meta_analysis_utils.r")
gene_name_mapping_file <- "../data/gene_mapping.txt.gz"

get_burden_and_SE <- function(
    file_paths,
    dt_gene_phenotype,
    filter_phenotype = NULL,
    case_control_threshold=100)
{
	files <- strsplit(file_paths, split=",")[[1]]
    folder <- gsub("(.*\\/)(.*)", "\\1", files[1])
    file <- gsub("(.*\\/)(.*)", "\\2", files[1])
    file_info_template <- extract_file_info(file)

    dt_genes <- fread(gene_name_mapping_file)
    names(dt_genes) <- c("Region", "geneID_version", "start", "stop",
        "chr", "gene_symbol")

    setkey(dt_genes, "Region")
    setkey(dt_gene_phenotype, "Region")
    dt_gene_phenotype <- merge(dt_gene_phenotype, dt_genes, all.x=TRUE)
    dt_gene_phenotype$max_MAF <- as.numeric(dt_gene_phenotype$max_MAF)
    setkeyv(dt_gene_phenotype, c("Region", "max_MAF", "Group"))

    dt_list <- list()
    for (f in files) {
        cat(f, "\n")
        # Get information about the files
        folder <- gsub("(.*\\/)(.*)", "\\1", f)
        file <- gsub("(.*\\/)(.*)", "\\2", f)
        file_info <- extract_file_info(file)
        checks(file_info, file_info_template)

        # Throw an error if the case count is too low
        if (file_info$binary) {
            file_info$n_cases <- as.integer(file_info$n_cases)
            file_info$n_controls <- as.integer(file_info$n_controls)
            if (file_info$n_cases < as.integer(case_control_threshold)) {
                next
            }
            if (file_info$n_controls < as.integer(case_control_threshold)) {
                next
            }
        } else if (file_info$n < as.integer(case_control_threshold)) {
            next
        }

        dt_list[[file]] <- merge(fread(f, key=c("Region", "max_MAF", "Group")), dt_gene_phenotype)
        dt_list[[file]]$dataset <- file_info$dataset
        dt_list[[file]]$ancestry <- file_info$ancestry
        dt_list[[file]]$binary <- file_info$binary
    }
    dt <- rbindlist(dt_list, use.names=TRUE, fill=TRUE)
    print(dt)
    dt[, `:=`(lower = BETA_Burden - SE_Burden, upper = BETA_Burden + SE_Burden)]
    return(dt)
}

# First, determine the collection of all traits that have a significant burden association
# (either 0.001, 1e-4, and any combination of pLoF and damaging missense). Then, loop over everything to 
# create the required plots.

# Oct 25 freeze on BMRC
meta_results <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100", pattern="tsv.gz", full.names=TRUE)

Groups_to_plot <- c(
    "damaging_missense_or_protein_altering",
    "pLoF",
    "pLoF;damaging_missense_or_protein_altering")
max_MAFs_to_plot <- c("0.001", "1e-04")
Pvalue_cutoff <- 2.5e-6

forest_plot <- list()
sexes <- list()
plot_count <- 0
for (f in meta_results) {
    cat(paste0(f, "\n"))
    dt_tmp <- fread(f) %>% filter(
        Group %in% Groups_to_plot,
        max_MAF %in% max_MAFs_to_plot,
        class == "Burden")
    forest_plot_tmp <- unique(dt_tmp %>%
        filter(Pvalue < Pvalue_cutoff) %>%
        select(Region, Group, max_MAF))
    phenotypeID <-  gsub(".*\\/([A-Za-z0-9]+)_.*", "\\1", f)
    sexes[[phenotypeID]] <- gsub(".*\\/([A-Za-z0-9]+)_([A-Z]+)_.*", "\\2", f)
    print(sexes)
    forest_plot_tmp$phenotypeID <- phenotypeID
    forest_plot[[phenotypeID]] <- data.table(forest_plot_tmp,
        key=c("Region", "Group", "max_MAF"))
    # Extract the BETA_burden and SE_burden for the inverse variance weighting for all of the
    # 'significant' associations.
    dt_tmp <- dt_tmp %>% filter(type == "Inverse variance weighted")
    dt_tmp <- data.table(dt_tmp %>% select(Region, Group, max_MAF, BETA_Burden, SE_Burden, Pvalue),
        key=c("Region", "Group", "max_MAF"))
    dt_tmp[, `:=`(lower = BETA_Burden - SE_Burden, upper = BETA_Burden + SE_Burden)]
    dt_tmp[, `:=`(ancestry = "Meta-analysis", dataset = "Meta-analysis")]
    forest_plot[[phenotypeID]] <- merge(forest_plot[[phenotypeID]], dt_tmp)
    plot_count <- plot_count + nrow(forest_plot[[phenotypeID]])
    cat(paste0("number of significant associations to plot: ", nrow(forest_plot[[phenotypeID]]), "\n"))
    cat(paste0("rolling total of plots: ",  plot_count, "\n"))
}

# Oct 25 freeze on BMRC
file_paths <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks", full.names=TRUE, recursive=TRUE, pattern="txt.gz")
dt_genes <- fread(gene_name_mapping_file)
names(dt_genes) <- c("Region", "geneID_version", "start", "stop",
        "chr", "gene_symbol")
setkey(dt_genes, "Region")
# Add the meta-analysis label
renaming_plot_biobank_list[["Meta-analysis"]] <- "Meta-analysis"
ancestry_levels <- c("AFR", "AMR", "EAS", "EUR", "SAS", "Meta-analysis")

ii <- 1
forest_plot_dt <- list()
pdf(file="forest_plots_gnomAD_colors_oct25_freeze.pdf", width=5.5, height=5)
for (name in names(forest_plot)) {
    print(name)
    files <- paste0(grep(paste0(".*cleaned/.*gene.*\\.", name, "\\.", ".*", "\\.", sexes[[name]], "\\..*gene.*"), file_paths, value=TRUE), collapse=",")
    dt <- get_burden_and_SE(files, forest_plot[[name]] %>% select(-c("BETA_Burden", "SE_Burden", "lower", "upper", "Pvalue")))
    dt_meta <- merge(forest_plot[[name]], dt_genes, all.x=TRUE)
    dt_meta[, binary:=dt$binary[1]]
    dt_meta[, Pvalue_Burden:=Pvalue]
    dt_meta[, Pvalue:=NULL]
    dt <- rbind(dt, dt_meta, fill=TRUE)
    relevel <- sort(names(renaming_plot_biobank_list), decreasing=TRUE)
    relevel <- c("Meta-analysis", relevel[-which(relevel == "Meta-analysis")])
    dt[, dataset := factor(dataset, levels = relevel, labels=renaming_plot_biobank_list[relevel])]
    setkeyv(dt, c("Region", "max_MAF", "Group"))
    plot_combinations <- unique(dt %>% select(Region, max_MAF, Group))
    plot_combinations <- data.table(plot_combinations, key = c("Region", "max_MAF", "Group"))

    if (nrow(plot_combinations) == 0) { next }
    for (i in 1:nrow(plot_combinations)) {
        dt_to_plot <- merge(dt, plot_combinations[i, ])
        dt_to_plot[, ancestry := factor(ancestry, levels = ancestry_levels)]
        p <- ggplot(dt_to_plot,
            aes(x = dataset, y = BETA_Burden, color = ancestry, group = paste(dataset, ancestry))) +
          geom_pointrange(position = position_dodge(width = 0.5), aes(ymin = lower, ymax = upper)) +
          geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
          labs(
            title = paste0(
                renaming_phenotype_list[name], "\n",
                dt_to_plot$gene_symbol[1], ": ",
                dt_to_plot$Region[1]),
            subtitle = paste0(paste0(renaming_plot_group_list[strsplit(dt_to_plot$Group[1], split=";")[[1]]], collapse=",\n"), ";\n",
                "max MAF = ", dt_to_plot$max_MAF[1]),
            x = "",
            y = TeX("$\\beta_{Burden}$")
            ) + coord_flip() + theme_bw() + labs(color = "Genetic\nancestry") +
            scale_color_manual(values = pop_colors)
        forest_plot_dt[[ii]] <- dt_to_plot
        ii <- ii+1
        print(p)
    }
}
dev.off()

forest_plot_dt <- rbindlist(forest_plot_dt, fill=TRUE)
fwrite(forest_plot_dt, sep='\t', quote=FALSE, file="forest_plot_data_oct25_freeze.tsv.gz")


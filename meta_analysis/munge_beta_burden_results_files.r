#!/bin/Rscript
library(data.table)
library(dplyr)
source("meta_analysis_utils.r")

# Renaming errors
renaming_list <- list(
	`damaging_missense_or_protein_altering` = c("damaging_missense"),
	`other_missense_or_protein_altering` = c("other_missense"),
	`pLoF;damaging_missense_or_protein_altering` = c("pLoF;damaging_missense"),
	`pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous` = c("pLoF;damaging_missense;other_missense;synonymous"),
	Cauchy = c(NA)
)

correct_names <- c("pLoF",
	"damaging_missense_or_protein_altering",
	"other_missense_or_protein_altering",
	"synonymous",
	"pLoF;damaging_missense_or_protein_altering",
	"pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous",
	"Cauchy")

# DEV: additional check to ensure that the correct names (above) are used in the variant level files
# DEV: check that either a single annotation file or 23 annotation files are provided
# DEV: write the results to disk

main <- function(args)
{
	gene_data_table <- extract_group_file_information(args$group_file_regexp)
	gene_files <- strsplit(args$gene_results_file_paths, split=",")[[1]]
	variant_files <- strsplit(args$variant_results_file_paths, split=",")[[1]]

    dt_list <- list()

    if (length(gene_files) != length(variant_files)) {
    	stop("Error: number of comma separated files doesn't match for gene and variant results")
    }

    for (i in 1:length(gene_files))
    {
        # Throw an error if the gene file is not there
        if (!file.exists(gene_files[i])) {
            stop(paste("File:", gene_files[i], "does not exist."))
        }
        # Throw an error if the variant file is not there
        if (!file.exists(variant_files[i])) {
            stop(paste("File:", variant_files[i], "does not exist."))
        }

        gene_file_info <- extract_file_info(gene_files[i])
        variant_file_info <- extract_file_info(variant_files[i])

        for (attribute in c(
            "dataset",
            "last_name",
            "analysis_name",
            "phenotype",
            "freeze_number",
            "sex", 
            "ancestry",
            "n_cases",
            "n_controls", 
            "software",
            "date",
            "binary")) {
        	if (gene_file_info[[attribute]] != variant_file_info[[attribute]]) {
        		stop(paste("Error:", attribute, "does not match between variant and gene results files"))
        	}
        }

        dt <- combine_gene_and_variant_information(variant_files[i], gene_data_table)
        combined_annotations <- extract_combined_annotations(gene_files[i])
        dt_burden <- extract_beta_burden_across_MAFs(dt, c(0.01, 0.001, 0.0001), combined_annotations)
        dt_burden <- dt_burden %>% 
            mutate(BETA_Burden_no_weighting = (numerator_common + numerator_rare) / (denominator_common + denominator_rare)) %>%
            select(Region, Group, max_MAF, BETA_Burden_no_weighting)

        dt_gene <- fread(gene_files[i])
        dt_burden <- data.table(dt_burden)
        setkeyv(dt_gene, c("Region", "Group", "max_MAF"))
        setkeyv(dt_burden, c("Region", "Group", "max_MAF"))
        dt_burden <- merge(dt_gene, dt_burden, all=TRUE)
        # Write the results
        # Remove redundant columns and relabel the beta_burden column
    }
}
	variant_files <- "data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.variant.20240110.txt.gz"
	gene_file <- "data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.gene.20240110.txt.gz"

	dt_burden <- extract_beta_burden_across_MAFs(dt, c(0.01, 0.001, 0.0001), combined_annotations)

# Then, extract the relevant variants from the variant results files
# Group in the appropriate way and apply the scaling to get the beta_burden to be used
# Merge this into the existing gene results, and apply the existing run_heterogeneity script


}

args <- list()
args$gene_results_file_paths <- "data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.AFR.175.6578.SAIGE.gene.20240110.cleaned.txt.gz,data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.gene.20240110.cleaned.txt.gz,data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.SAS.610.6233.SAIGE.gene.20240110.cleaned.txt.gz"
args$variant_results_file_paths <- "data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.AFR.175.6578.SAIGE.variant.20240110.txt.gz,data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.EUR.19915.382460.SAIGE.variant.20240110.txt.gz,data/meta_analysis/gcloud/uk-biobank.palmer.PRELIMINARY.Coronary_artery_disease.JULY23Freeze.ALL.SAS.610.6233.SAIGE.variant.20240110.txt.gz"
args$group_file_regexp <- "ukb_wes_450k.july.qced.brava_common_rare.v7.chr.*.saige.txt.gz"

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--gene_results_file_paths", default=NULL, required=TRUE,
    help="Comma separated SAIGE gene results file paths for files for meta-analysis")
parser$add_argument("--variant_results_file_paths", default=NULL, required=TRUE,
    help="Comma separated SAIGE variant results file paths for files for meta-analysis")
parser$add_argument("--group_file_regexp", default=NULL, required=TRUE,
    help="Regular expression for the group files to combine for this biobank")
args <- parser$parse_args()

main(args)

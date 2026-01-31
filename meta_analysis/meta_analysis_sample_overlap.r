#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

# Pass a comma separated list of filenames
tests <- c("Burden", "SKAT", "SKAT-O")
source("meta_analysis_utils.r")

# Extract the files to read in
# Pass the name of the trait to be used in the output
# Let the case control threshold be an option
# Let the N threshold be an option

main <- function(args)
{
    files <- strsplit(args$file_paths, split=",")[[1]]
    case_control_threshold <- as.integer(args$case_control_threshold)
    dt_list <- list()
    folder <- gsub("(.*\\/)(.*)", "\\1", files[1])
    file <- gsub("(.*\\/)(.*)", "\\2", files[1])
    file_info_template <- extract_file_info(file)
    cat(paste(paste(names(args), args, sep=": "), collapse="\n"), "\n")

    for (f in files)
    {
        # Throw an error if the file is not there
        if (!file.exists(f)) {
            stop(paste("File:", file, "does not exist."))
        }

        # Get information about the file
        folder <- gsub("(.*\\/)(.*)", "\\1", f)
        file <- gsub("(.*\\/)(.*)", "\\2", f)
        file_info <- extract_file_info(file)
        checks(file_info, file_info_template, args$no_sex_check)

        # Throw an error if the case count is too low
        if (file_info$binary) {
            file_info$n_cases <- as.integer(file_info$n_cases)
            file_info$n_controls <- as.integer(file_info$n_controls)
            if (file_info$n_cases < case_control_threshold) {
                next
            }
            if (file_info$n_controls < case_control_threshold) {
                next
            }
        } else {
            file_info$n <- as.integer(file_info$n)
            if (file_info$n < case_control_threshold) {
                next
            }
        }

        dt_list[[file]] <- fread(paste0(folder, file))
        dt_list[[file]]$dataset <- file_info$dataset
        dt_list[[file]]$ancestry <- file_info$ancestry

        # Compare N to actual N and throw an error if it doesn't match
        dt_tmp <- add_N_using_filename(file_info, dt_list[[file]])
         if (!args$no_Neff) {
            dt_tmp <- add_N_using_Neff_weights_file(file_info, dt_tmp,
            Neff_weights_file=args$Neff_weights_file)
        }
        dt_list[[file]] <- dt_tmp %>% filter(Group != "Cauchy")
    }

    dt <- rbindlist(dt_list, use.names=TRUE)
    dt_cor <- determine_null_correlation(dt, binary=file_info$binary) # Note that this is just using the Burden p-values
    # Merge with Neff information
    dt <- dt %>% filter((!is.na(Pvalue)) & (!is.na(Pvalue_SKAT)) & (!is.na(Pvalue_Burden)))

    # Nudge P-values that are close to exactly 1 (for Stouffer).
    dt <- dt %>% mutate(
        Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue),
        Pvalue_SKAT = ifelse(Pvalue_SKAT > 0.99, 0.99, Pvalue_SKAT),
        Pvalue_Burden = ifelse(Pvalue_Burden > 0.99, 0.99, Pvalue_Burden)
    )
    setkey(dt, "dataset")

    dt_n_eff <- unique(dt %>% select(Region, dataset, ancestry, N_eff))
    setkeyv(dt_n_eff, c("Region", "dataset", "ancestry", "N_eff"))
    cauchy_only <- ifelse(grepl("\\.extra_cauchy\\.", files[1]), TRUE, FALSE)

    if (!cauchy_only) {
        dt_meta <- list()
        for (test in tests)
        {
            cat(paste0(test, "...\n"))
            dt_meta[[test]] <- list()
            Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))
            # Lets determine what the additional cross term column will be - this will be
            # what is summed over in the denominator in Stouffer
            dt <- dt %>%
                group_by(Region, Group, max_MAF, ancestry) %>%
                mutate(
                    n_rows = n(),
                    # Calculate cross_terms using expand.grid within each group
                    cross_terms = {
                        dataset_pairs <- cbind(expand.grid(dataset1 = dataset, dataset2 = dataset), ancestry=ancestry)
                        cross_sum <- dt_cor[.(dataset_pairs)]$summation
                        cross_sum_total <- colSums(matrix(cross_sum, nrow = n_rows), na.rm = TRUE)
                        # If cross_sum_total is zero, use N_eff; otherwise, use cross_sum_total
                        ifelse(cross_sum_total == 0, N_eff, cross_sum_total)
                        }
                ) %>% ungroup()

            # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
            dt_meta[[test]][["Stouffer"]] <- run_stouffer_overlap(dt %>% group_by(Region, Group, max_MAF),
                "N_eff", "Stat", "cross_terms", Pvalue_col, "Pvalue",
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(type="Stouffer")
            dt_meta[[test]][["Stouffer"]] <- data.table(
                dt_meta[[test]][["Stouffer"]], key = c("Region", "Group", "max_MAF"))

            dt_meta[[test]] <- rbindlist(dt_meta[[test]], use.names=TRUE, fill=TRUE) %>% mutate(class=test)
        }

        dt_meta <- rbindlist(dt_meta, fill=TRUE)

        # This is evaluating Cauchy ourselves, for everything
        dt_cauchy <- list()
        for (test in tests) {
            cat(paste0(test, " Cauchy combination..."))
            Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))        
            dt_cauchy[[test]] <- run_cauchy(dt %>% group_by(Region, dataset, ancestry),
                "N_eff", "Stat", Pvalue_col, "Pvalue") %>% mutate(type="Cauchy")
            dt_cauchy[[test]] <- merge(dt_cauchy[[test]], dt_n_eff) %>% mutate(Group = "Cauchy")
        }
        cat("\n")
    } else {
        dt_cauchy <- list()
        for (test in tests) {
            cat(paste0(test, " Cauchy combination..."))
            Pvalue_col <- ifelse(test == "SKAT-O", "Pvalue", paste0("Pvalue_", test))  
            dt_cauchy[[test]] <- dt %>% 
                select(Region, Group, dataset, ancestry, .data[[Pvalue_col]], N_eff) %>%
                rename(Pvalue = .data[[Pvalue_col]])
            dt_cauchy[[test]] <- data.table(dt_cauchy[[test]], key=c("Region", "dataset", "ancestry", "N_eff"))
            dt_cauchy[[test]] <- merge(dt_cauchy[[test]], dt_n_eff)
        }
    }

    dt_meta_cauchy <- list()
    for (test in tests)
    {
        dt_cauchy[[test]] <- dt_cauchy[[test]] %>%
            group_by(Region, Group, ancestry) %>%
            mutate(
                n_rows = n(),
                # Calculate cross_terms using expand.grid within each group
                cross_terms = {
                    dataset_pairs <- cbind(expand.grid(dataset1 = dataset, dataset2 = dataset), ancestry=ancestry)
                    cross_sum <- dt_cor[.(dataset_pairs)]$summation
                    cross_sum_total <- colSums(matrix(cross_sum, nrow = n_rows), na.rm = TRUE)
                    # If cross_sum_total is zero, use N_eff; otherwise, use cross_sum_total
                    ifelse(cross_sum_total == 0, N_eff, cross_sum_total)
                    }
            ) %>% ungroup()
        dt_meta_cauchy[[test]] <- list()
        # Stouffer's Z - Make sure P-values match, Stat = weighted_Z_Burden_Stouffer
        dt_meta_cauchy[[test]][["Stouffer"]] <- run_stouffer_overlap(
            dt_cauchy[[test]] %>% group_by(Region, Group),
            "N_eff", "Stat", "cross_terms", "Pvalue", "Pvalue") %>% mutate(type="Stouffer")
        dt_meta_cauchy[[test]] <- rbindlist(dt_meta_cauchy[[test]], use.names=TRUE) %>% mutate(class=test)
    }
    dt_meta_cauchy <- rbindlist(dt_meta_cauchy)

    if (!cauchy_only)
    {
        cat("writing the results...\n")
        dt_meta <- rbind(dt_meta, dt_meta_cauchy %>% mutate(Group="Cauchy", max_MAF="Cauchy"), fill=TRUE)
        fwrite(dt_meta, file=ifelse(grepl(".tsv.gz$", args$out), args$out, paste0(args$out, ".tsv.gz")), sep='\t')
    } else {
        cat("writing the results...\n")
        fwrite(dt_meta_cauchy, file=ifelse(grepl(".tsv.gz$", args$out), args$out, paste0(args$out, ".tsv.gz")), sep='\t')
    }
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--file_paths", default=NULL, required=TRUE,
    help="Comma separated file paths for files for meta-analysis")
parser$add_argument("--case_control_threshold", default=100, required=FALSE,
    help=paste0("Case/control count threshold for inclusion of the data into the meta-analysis [default=100]"))
parser$add_argument("--out", default="meta_analysis", required=FALSE,
    help="Output filepath")
parser$add_argument("--no_sex_check", default=FALSE, action='store_true',
    help="Perform sex check of samples used for the trait when running meta-analysis?")
parser$add_argument("--Neff_weights_file",
    default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_may25.tsv.gz",
    help="File to pass effective sample sizes")
parser$add_argument("--no_Neff", default=FALSE, action='store_true',
    help="run using Neff estimated from case-control counts rather than the GRM")
args <- parser$parse_args()

main(args)

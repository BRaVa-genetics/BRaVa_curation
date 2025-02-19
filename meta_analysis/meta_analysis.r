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

        print(names(dt_list[[file]]))

        # Compare N to actual N and throw an error if it doesn't match
        dt_tmp <- add_N_using_filename(file_info, dt_list[[file]])
        dt_tmp <- add_N_using_Neff_weights_file(file_info, dt_tmp,
            Neff_weights_file=args$Neff_weights_file)
        dt_list[[file]] <- dt_tmp %>% filter(Group != "Cauchy")
    }

    print(dt_list)
    for (i in 1:length(dt_list)) {
        print(names(dt_list)[i])
        print(names(dt_list[[i]]))
        print(length(names(dt_list[[i]])))
    }
    dt <- rbindlist(dt_list, use.names=TRUE, fill=TRUE)
    dt <- dt %>% filter((!is.na(Pvalue)) & (!is.na(Pvalue_SKAT)) & (!is.na(Pvalue_Burden)))

    # Nudge P-values that are close to exactly 1 (for Stouffer).
    dt <- dt %>% mutate(
        Pvalue = ifelse(Pvalue > 0.99, 0.99, Pvalue),
        Pvalue_SKAT = ifelse(Pvalue_SKAT > 0.99, 0.99, Pvalue_SKAT),
        Pvalue_Burden = ifelse(Pvalue_Burden > 0.99, 0.99, Pvalue_Burden)
    )

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
            
            # Weighted Fisher's meta-analysis of p-values
            dt_meta[[test]][["weighted Fisher"]] <- run_weighted_fisher(
                dt %>% group_by(Region, Group, max_MAF),
                "N_eff", Pvalue_col, "Pvalue",
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(Stat = NA, type="Weighted Fisher")

            # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
            dt_meta[[test]][["Stouffer"]] <- run_stouffer(dt %>% group_by(Region, Group, max_MAF),
                "N_eff", "Stat", Pvalue_col, "Pvalue",
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(type="Stouffer")
            dt_meta[[test]][["Stouffer"]] <- data.table(
                dt_meta[[test]][["Stouffer"]], key = c("Region", "Group", "max_MAF"))

            if (test == "Burden") {
                # And also run the inverse-variance weighted meta-analysis
                dt_meta[[test]][["inverse_variance_weighted"]] <- run_inv_var(
                    dt %>% group_by(Region, Group, max_MAF), "BETA_Burden", "SE_Burden",
                    "BETA_Burden", "SE_Burden", "Pvalue") %>%
                mutate(type="Inverse variance weighted")
                dt_meta[[test]][["inverse_variance_weighted"]] <- data.table(
                    dt_meta[[test]][["inverse_variance_weighted"]], key = c("Region", "Group", "max_MAF"))

                # And also evaluate the heterogeneity P-values for both inverse variance weighted
                # and Stouffers
                dt_meta[[test]][["Stouffer"]] <- merge(dt_meta[[test]][["Stouffer"]],
                    data.table(run_heterogeneity_test(
                        weights(dt, FALSE, se_name="SE_Burden", n_eff_name="N_eff") %>% 
                            group_by(Region, Group, max_MAF),
                        input_beta="BETA_Burden", output_meta_beta="BETA_Burden"),
                    key=c("Region", "Group", "max_MAF"))
                )
                dt_meta[[test]][["inverse_variance_weighted"]] <- merge(dt_meta[[test]][["inverse_variance_weighted"]],
                    data.table(run_heterogeneity_test(
                        weights(dt, TRUE, se_name="SE_Burden", n_eff_name="N_eff") %>% 
                            group_by(Region, Group, max_MAF),
                        input_beta="BETA_Burden", output_meta_beta="BETA_meta"),
                    key=c("Region", "Group", "max_MAF"))
                )
            }
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
        dt_meta_cauchy[[test]] <- list()
        # Weighted Fisher's meta-analysis of p-values
        dt_meta_cauchy[[test]][["weighted Fisher"]] <- run_weighted_fisher(
            dt_cauchy[[test]] %>% group_by(Region, Group),
            "N_eff", "Pvalue", "Pvalue") %>% mutate(Stat = NA, type="Weighted Fisher")
        # Stouffer's Z - Make sure P-values match, Stat = weighted_Z_Burden_Stouffer
        dt_meta_cauchy[[test]][["Stouffer"]] <- run_stouffer(dt_cauchy[[test]] %>% group_by(Region, Group),
            "N_eff", "Stat", "Pvalue", "Pvalue") %>% mutate(type="Stouffer")
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
    default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights.tsv.gz",
    help="File to pass effective sample sizes")
args <- parser$parse_args()

main(args)

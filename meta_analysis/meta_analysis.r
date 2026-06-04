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

# estimate_var <- function(CAF, n_rare, n_ultra_rare)
# {
#     lambda <- 2 * CAF
#     if (lambda < 0.05) {
#         var <- 2 * CAF * (1 - CAF)
#     } else if ((lambda < 0.2) & (n_rare > 0)) {
#         var <- 2 * (CAF - CAF^2 / m)
#     } else {
#         frac_rare  <- Number_rare / (Number_rare + Number_ultra_rare)
#         CAF_rare   <- CAF * frac_rare
#         CAF_ultra  <- CAF * (1 - frac_rare)

#         if (n_rare > 0) {
#             var_rare <- 2 * (CAF_rare - CAF_rare^2 / Number_rare)
#         } else {
#             var_rare <- 0
#         }

#         lambda_u <- 2 * CAF_ultra
#         c <- 1 - exp(-lambda_u)
#         var_ultra <- c * (1 - c)
#         var <- var_rare + var_ultra
#     }
#     return(var)
# }

estimate_var <- function(CAF, n_rare, n_ultra_rare) {

  lambda <- 2 * CAF
  var <- numeric(length(CAF))

  # Case 1: very small lambda
  idx1 <- lambda < 0.05
  var[idx1] <- 2 * CAF[idx1] * (1 - CAF[idx1])

  # Case 2: intermediate lambda with rare variants
  idx2 <- (lambda >= 0.05) & (lambda < 0.2) & (n_rare > 0)
  var[idx2] <- 2 * (CAF[idx2] - CAF[idx2]^2 / n_rare[idx2])

  # Case 3: everything else (split rare / ultra-rare)
  idx3 <- !(idx1 | idx2)

  if (any(idx3)) {
    frac_rare <- n_rare[idx3] / (n_rare[idx3] + n_ultra_rare[idx3])
    frac_rare[is.na(frac_rare)] <- 0

    CAF_rare  <- CAF[idx3] * frac_rare
    CAF_ultra <- CAF[idx3] * (1 - frac_rare)

    var_rare <- numeric(sum(idx3))
    has_rare <- n_rare[idx3] > 0
    var_rare[has_rare] <-
      2 * (CAF_rare[has_rare] -
           CAF_rare[has_rare]^2 / n_rare[idx3][has_rare])

    lambda_u <- 2 * CAF_ultra
    c <- 1 - exp(-lambda_u)
    var_ultra <- c * (1 - c)

    var[idx3] <- var_rare + var_ultra
  }

  return(pmax(var, 0))
}

main <- function(args)
{
    files <- strsplit(args$file_paths, split=",")[[1]]
    case_control_threshold <- as.integer(args$case_control_threshold)
    dt_list <- list()
    folder <- gsub("(.*\\/)(.*)", "\\1", files[1])
    file <- gsub("(.*\\/)(.*)", "\\2", files[1])
    file_info_template <- extract_file_info(file)
    cat(paste(paste(names(args), args, sep=": "), collapse="\n"), "\n")
    if (args$estimate_CAF) {
        ref <- fread(args$reference_CAF)
        setkeyv(ref, c("Region", "ancestry", "max_MAF", "Group", "dataset"))
        ref_max <- ref[, .SD[n_ref == max(as.integer(n_ref))], by = ancestry]
        ref_max <- ref_max %>% select(-dataset)
        setkeyv(ref_max, c("Region", "ancestry", "max_MAF", "Group"))
    }

    for (f in files)
    {
        print(f)
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
        dt_tmp <- dt_tmp %>% filter(Group != "Cauchy")
        setDT(dt_tmp)

        if (args$estimate_CAF) {
            if (file_info$binary) {
                if ("MAC_case" %in% names(dt_tmp)) {
                    # No need to impute - determine estimate directly
                    cat(paste0(file_info$dataset, ", ", file_info$ancestry, ": No need to impute - determine estimate directly\n"))
                    dt_tmp[, CAF := MAC / (2*(file_info$n_cases+file_info$n_controls))]
                } else {
                    cat(paste0(file_info$dataset, ", ", file_info$ancestry, ": Impute!\n"))
                    setindexv(dt_tmp, c("Region", "max_MAF", "Group"))
                    if (ref[ancestry == file_info$ancestry & dataset == file_info$dataset, .N] > 0) {
                        cat(paste0("using: ", file_info$dataset, ": ", file_info$ancestry, " for imputation\n"))
                        dt_tmp <- ref[ancestry == file_info$ancestry & dataset == file_info$dataset, ][dt_tmp, on = .(Region, max_MAF, Group)]
                        dt_tmp[is.na(CAF), CAF := 1 / (2*(file_info$n_cases+file_info$n_controls))]
                        if ("i.Number_rare" %in% names(dt_tmp)) {
                            dt_tmp[, `:=`(
                              Number_rare  = i.Number_rare,
                              Number_ultra_rare = i.Number_ultra_rare
                            )]
                            dt_tmp[, c("i.Number_rare", "i.Number_ultra_rare") := NULL]
                        }
                        dt_tmp[is.na(Number_ultra_rare), Number_ultra_rare := 1]
                        dt_tmp[is.na(Number_rare), Number_rare := 0]
                        dt_tmp[, `:=`(
                          dataset  = i.dataset,
                          ancestry = i.ancestry
                        )]
                        dt_tmp[, c("i.dataset", "i.ancestry") := NULL]
                    } else {
                        cat("No exact match, using largest cohort in the references\n")
                        dt_tmp <- ref_max[ancestry == file_info$ancestry, ][dt_tmp, on = .(Region, max_MAF, Group)]
                        dt_tmp[is.na(CAF), CAF := 1 / (2*(file_info$n_cases+file_info$n_controls))]
                        if ("i.Number_rare" %in% names(dt_tmp)) {
                            dt_tmp[, `:=`(
                              Number_rare  = i.Number_rare,
                              Number_ultra_rare = i.Number_ultra_rare
                            )]
                            dt_tmp[, c("i.Number_rare", "i.Number_ultra_rare") := NULL]
                        }
                        dt_tmp[is.na(Number_ultra_rare), Number_ultra_rare := 1]
                        dt_tmp[is.na(Number_rare), Number_rare := 0]
                        dt_tmp[, ancestry := i.ancestry]
                        dt_tmp[, `i.ancestry` := NULL]
                    }
                }
            } else {
                if ("MAC" %in% names(dt_tmp)) {
                    # No need to impute - determine estimate directly
                    cat(paste0(file_info$dataset, ", ", file_info$ancestry, ": No need to impute - determine estimate directly\n"))
                    dt_tmp[, CAF := MAC / (2*(file_info$n))]
                } else {
                    cat(paste0(file_info$dataset, ", ", file_info$ancestry, ": Impute!\n"))
                    setindexv(dt_tmp, c("Region", "max_MAF", "Group"))
                    if (ref[ancestry == file_info$ancestry & dataset == file_info$dataset, .N] > 0) {
                        cat(paste0("using: ", file_info$dataset, ": ", file_info$ancestry, " for imputation\n"))
                        dt_tmp <- ref[ancestry == file_info$ancestry & dataset == file_info$dataset, ][dt_tmp, on = .(Region, max_MAF, Group)]
                        dt_tmp[is.na(CAF), CAF := 1 / (2*file_info$n)]
                        if ("i.Number_rare" %in% names(dt_tmp)) {
                            dt_tmp[, `:=`(
                              Number_rare  = i.Number_rare,
                              Number_ultra_rare = i.Number_ultra_rare
                            )]
                            dt_tmp[, c("i.Number_rare", "i.Number_ultra_rare") := NULL]
                        }
                        dt_tmp[is.na(Number_ultra_rare), Number_ultra_rare := 1]
                        dt_tmp[is.na(Number_rare), Number_rare := 0]
                        dt_tmp[, `:=`(
                          dataset  = i.dataset,
                          ancestry = i.ancestry
                        )]
                        dt_tmp[, c("i.dataset", "i.ancestry") := NULL]
                    } else {
                        cat("No exact match, using largest cohort in the references\n")
                        dt_tmp <- ref_max[ancestry == file_info$ancestry, ][dt_tmp, on = .(Region, max_MAF, Group)]
                        dt_tmp[is.na(CAF), CAF := 1 / (2*file_info$n)]
                        if ("i.Number_rare" %in% names(dt_tmp)) {
                            dt_tmp[, `:=`(
                              Number_rare  = i.Number_rare,
                              Number_ultra_rare = i.Number_ultra_rare
                            )]
                            dt_tmp[, c("i.Number_rare", "i.Number_ultra_rare") := NULL]
                        }
                        dt_tmp[is.na(Number_ultra_rare), Number_ultra_rare := 1]
                        dt_tmp[is.na(Number_rare), Number_rare := 0]
                        dt_tmp[, ancestry := i.ancestry]
                        dt_tmp[, `i.ancestry` := NULL]
                    }
                }
            }
        }

        dt_list[[file]] <- dt_tmp
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

            if (args$estimate_CAF) {
                dt_to_test <- dt %>% 
                    mutate(var=estimate_var(CAF, Number_rare, Number_ultra_rare)) %>% 
                    mutate(Neff_caf=N_eff*var) %>% 
                    group_by(Region, Group, max_MAF)
            } else {
                dt_to_test <- dt %>% group_by(Region, Group, max_MAF)
            }

            # Weighted Fisher's meta-analysis of p-values
            dt_meta[[test]][["weighted Fisher"]] <- run_weighted_fisher(
                dt_to_test, "N_eff", Pvalue_col, "Pvalue",
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(Stat = NA, type="Weighted Fisher")

            # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
            dt_meta[[test]][["Stouffer"]] <- run_stouffer(
                dt_to_test, "N_eff", "Stat", Pvalue_col, "Pvalue",
                two_tail = ifelse(test == "Burden", TRUE, FALSE),
                input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
            mutate(type="Stouffer")
            dt_meta[[test]][["Stouffer"]] <- data.table(
                dt_meta[[test]][["Stouffer"]], key = c("Region", "Group", "max_MAF"))

            if (args$estimate_CAF) {
                # Weighted Fisher's meta-analysis of p-values
                dt_meta[[test]][["weighted Fisher CAF"]] <- run_weighted_fisher(
                    dt_to_test, "Neff_caf", Pvalue_col, "Pvalue",
                    two_tail = ifelse(test == "Burden", TRUE, FALSE),
                    input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
                mutate(Stat = NA, type="Weighted Fisher CAF")

                # Stouffer's Z - Make sure P-values match, Stat= weighted_Z_Burden_Stouffer
                dt_meta[[test]][["Stouffer CAF"]] <- run_stouffer(
                    dt_to_test, "Neff_caf", "Stat", Pvalue_col, "Pvalue",
                    two_tail = ifelse(test == "Burden", TRUE, FALSE),
                    input_beta = ifelse(test == "Burden", "BETA_Burden", NULL)) %>% 
                mutate(type="Stouffer CAF")
                dt_meta[[test]][["Stouffer CAF"]] <- data.table(
                    dt_meta[[test]][["Stouffer CAF"]], key = c("Region", "Group", "max_MAF"))
            }

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
    default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/Neff/Neff_weights_may25.tsv.gz",
    help="File to pass effective sample sizes")
parser$add_argument("--no_Neff", default=FALSE, action='store_true',
    help="run using Neff estimated from case-control counts rather than the GRM")
parser$add_argument("--estimate_CAF", default=FALSE, action='store_true',
    help="Estimate CAF and use in additional weighted Stouffer analysis")
parser$add_argument("--reference_CAF", default="/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/full_reference_cafs.tsv.gz",
    help="File containing CAF references to use for imputation for CAF weighted Stouffer", required=FALSE)
args <- parser$parse_args()

main(args)

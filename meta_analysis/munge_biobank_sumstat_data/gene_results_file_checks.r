#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

source("../meta_analysis_utils.r")

main <- function(args)
{
    files <- strsplit(args$file_paths, split=",")[[1]]
    for (f in files)
    {
        cat(paste0("checking file...", f, "...\n"))
        
        # Throw an error if the file is not there
        if (!file.exists(f) & !args$just_filename) {
            stop(paste("File:", file, "does not exist."))
        }

        # Get information about the file
        folder <- gsub("(.*\\/)(.*)", "\\1", f)
        file <- gsub("(.*\\/)(.*)", "\\2", f)
        file_info <- extract_file_info(file)

        cat(paste0("\nfile information extracted from file naming according to specified convention:\n",
            "continuous traits: [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz\n",
            "binary traits: [dataset].[last_name].[analysis_name].[phenotype].[freeze_number].[sex].[ancestry].[n_cases].[n_controls].[SAIGE].{gene,variant}.[YYYYMMDD].txt.gz\n\n"
        ))
        cat(paste(names(file_info), file_info, sep = ": "), sep = "\n")

        cat("\nchecking file naming...\n")
        # Use this to perform series of checks
        if (!(file_info$dataset %in% names(file_check_information$dataset))) {
            cat(paste0("accepted dataset names:\n", paste(names(file_check_information$dataset), collapse="\n"), "\n"))
            stop(paste0(file_info$dataset, " is not in the accepted collection of dataset names"))
        }

        if (!(file_info$phenotype %in% names(file_check_information$phenotype))) {
            cat(paste0("accepted (unique) phenotype abbreviations:\n", paste(file_check_information$phenotype, collapse="\n"), "\n"))
            stop(paste0(file_info$phenotype, " is not in the accepted collection of dataset names"))
        }

        if (!(file_info$ancestry %in% names(file_check_information$ancestry))) {
            cat(paste0("accepted population labels:\n", paste(file_check_information$ancestry, collapse="\n"), "\n"))
            stop(paste0(file_info$ancestry, " is not in the accepted collection of population labels"))
        }

        if (!(file_info$sex %in% names(file_check_information$sex))) {
            cat(paste0("accepted sex labels:\n", paste(file_check_information$sex, collapse="\n"), "\n"))
            stop(paste0(file_info$sex, " is not in the accepted collection of sex labels"))
        }

        if (file_info$type != args$type) {
            cat(paste0("This checking script is expecting ", args$type, " level results\n"))
            stop(paste0(file_info$type, " is not '", args$type, "'"))
        }

        cat("naming convention checking complete!\n")

        if (!file_info$gz) {
            cat("\nThe file does not appear to be compressed\n",
                "consider compressing before uploading to the cloud\n")
        }

        if (args$just_filename) {
             cat("just file naming...\n")
             break
        }

        dt_header <- fread(cmd = ifelse(file_info$gz,
            paste0("gzcat ", folder, file),
            paste0("cat ", folder, file)), nrows=0)

        results_file_columns <- names(dt_header)

        if (length(results_file_columns) > length(unique(results_file_columns))) {
            stop(paste0(
                "multiple columns with the same name are present within the ",
                args$type, " results file"))
        }

        if (args$type == "variant")
        {
            # Now, ensure that all of the required columns are present
            if (all(default_variant_result_columns[[ifelse(file_info$binary, "binary", "continuous")]] %in% results_file_columns))
            {
                cat("\nall required gene results columns are present\n")
            } else {
                cat("\na collection of columns are missing, perhaps there is a naming issue?\n")
                cat(paste0("columns missing:\n",
                    setdiff(default_variant_result_columns[[ifelse(file_info$binary, "binary", "continuous")]], results_file_columns), collapse="\n"))
                stop("required columns missing from the results file")
            }
        }

        if (args$type == "gene") 
        {
            # Now, ensure that all of the required columns are present
            if (all(default_gene_result_columns$minimal %in% results_file_columns))
            {
                cat("\nall required gene results columns are present\n")
                if (file_info$binary) {
                    if (all(default_gene_result_columns$binary %in% results_file_columns)) {
                        cat("allele counts are also present\n")
                    }
                }
                if (!file_info$binary) {
                    if (all(default_gene_result_columns$continuous %in% results_file_columns)) {
                        cat("allele counts are also present\n")
                    }
                }

            } else {
                cat("\na collection of columns are missing, perhaps there is a naming issue?\n")
                cat(paste0("columns missing:\n",
                    setdiff(minimal_gene_result_columns, results_file_columns), collapse="\n"))
                stop("required columns missing from the results file")
            }

            # Determine that all requested annotations are present
            dt <- fread(cmd = ifelse(file_info$gz,
                paste0("gzcat ", folder, file),
                paste0("cat ", folder, file)))

            annotations <- unique(dt$Group)
            if (all(correct_names %in% annotations)) {
                cat("\nall requested annotations are present\n")
            } else {
                cat("a subset of the requested annotations are not present\n")
                cat("\nrequested annotations:\n", paste(minimal_group_list, collapse="\n "),"\n")
                cat("\nannotations present:\n", paste(annotations, collapse="\n "), "\n")
                cat("\nannotations missing:\n",
                    paste0(setdiff(minimal_group_list, annotations), collapse="\n "), "\n")
                cat("\nthis could be due to a naming difference: please edit the naming in your results file\n")
                stop("a subset of the requested annotations are not present\n")
            }

            max_MAFs <- unique(dt$max_MAF)

            if (all(max_MAFs %in% requested_max_MAFs)) {
                cat("\nall requested max_MAF thresholds are present")
            } else {
                if (length(max_MAFs > 0) & 
                    (max(max_MAFs, na.rm=TRUE) < max(requested_max_MAFs))
                ) {
                    stop("Not all max MAF thresholds are present, a higher max_MAF threshold was not run")
                } else {
                    cat(paste0("\nWarning: not all max_MAF thresholds were run; ",
                        "this is likely because this subset of the dataset you are looking at is quite small\n"))
                }
            }

            # Determine that the format of the gene names is correct
            if (all(grepl("ENSG[0-9]+", dt$Region))) {
                cat("\nall Region names contain ensembl gene IDs\n")
            } else {
                stop("please use ensembl gene IDs only for the region labels")
            }
        }
        
        cat(paste0("file: ", f,
            " appears to be correctly named and contains the information we need!\n",
            "Thank you Prof/Dr/Ms/Mr/Miss ", file_info$last_name, "!\n"))
    }
    cat("all files passed!\n")
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--file_paths", default=NULL, required=TRUE,
    help="Comma separated file paths for gene-based results files for meta-analysis")
parser$add_argument("--just_filename", default=FALSE, action="store_true",
    help="Comma separated file paths for gene-based results files for meta-analysis")
parser$add_argument("--type", default="gene", required=FALSE)
args <- parser$parse_args()

main(args)

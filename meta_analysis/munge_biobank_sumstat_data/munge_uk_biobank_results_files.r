#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

# For naming files
biobank <- "uk-biobank"
last_name <- "palmer"
analysis_name <- "pilot"
freeze_number <- "JULY23Freeze"
date <- "20240110"
method <- "SAIGE"

dx_data_dir <- "brava/outputs/step2/sept2023/"
data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/raw")
out_data_dir <- paste0("~/Repositories/BRaVa_curation/data/meta_analysis/gcloud/", biobank, "/cleaned")
system(paste0("mkdir -p ", data_dir, "/gene/sept2023/combined"))
system(paste0("mkdir -p ", data_dir, "/variant/sept2023/combined"))
system(paste0("mkdir -p ", out_data_dir, "/gene"))
system(paste0("mkdir -p ", out_data_dir, "/variant"))
download <- TRUE

# DEV: remember to apply the checks to all of the munge biobank X scripts
# DEV: we may require an additional loop restricting to sex specific analysis in the future.
# Note, this is for round 1 of the uk-biobank runs.
# We will rerun the analysis, restricting to Barney's QCed samples and variants 
# and superpopulation labellings as well, and compare results.
download_cleaned <- TRUE

if (!download_cleaned) {
    if (download) {
    	# Replace with dx download, and we now need to combine across chromosomes
    	# Also, we should flag if a chromosome is missing!
    	system(paste0(
    		"dx download -r ", dx_data_dir, "/ -o ", data_dir , "/gene/")
    	)
    	system(paste0("mv ", data_dir, "/gene/sept2023/*singleAssoc* ", data_dir, "/variant/sept2023/"))
    }

    raw_output_folder <- paste0(data_dir, "/gene/sept2023")
    gene_files <- grep(".gz$", dir(raw_output_folder, full.names=TRUE), value=TRUE)
    phenotypes <- unique(gsub(".*chr[0-9X]+_(split_[0-9]+_)*(.*)_(AFR|AMR|EAS|EUR|SAS).*", "\\2", gene_files))

    # Extract additional information for file naming
    # (we also require the variant files to obtain the counts for the file names, 
    # extract a single chromosome file for each of the variant output files and determine the number
    # of case and controls)

    raw_output_variant_folder <- paste0(data_dir, "/variant/sept2023")
    variant_files <- grep(".gz$", dir(raw_output_variant_folder, full.names=TRUE), value=TRUE)
    phenotype_info <- list()

    files <- list(gene=gene_files, variant=variant_files)

    for (phenotype in phenotypes) {
        cat(paste0(phenotype, ": "))
        subfiles <- grep(phenotype, variant_files, value=TRUE)
        pops <- unique(gsub(".*chr[0-9X]+_(split_[0-9]+_)*(.*)_(AFR|AMR|EAS|EUR|SAS).*", "\\3", subfiles))
        phenotype_info[[phenotype]] <- list()
        for (pop in pops) {
            cat(paste0(pop, ".."))
            subsubfiles <- grep(paste0("_", pop, "[_]?[\\.]?"), subfiles, value=TRUE)
            sex <- ifelse(grepl("_(F|M).txt", subsubfiles[1]), gsub(".*_(F|M).txt.*", "\\1", subsubfiles), "ALL")
            # Ensure that all chromosomes are present
            to_match <- paste0("chr",c(seq(1,22), "X"))
            first_match <- subsubfiles[which(to_match %in% unique(gsub(".*(chr[0-9X]+).*", "\\1", subsubfiles)))[1]]
            header <- fread(first_match, nrows=1)
            if ("N_case" %in% names(header)) {
                binary <- TRUE
                n_cases <- header$N_case[1]
                n_controls <- header$N_ctrl[1]
                phenotype_info[[phenotype]][[pop]] <- list(
                    binary=binary, N_case=n_cases, N_ctrl=n_controls, sex=sex)
            } else {
                binary <- FALSE
                n <- header$N
                phenotype_info[[phenotype]][[pop]] <- list(binary=binary, N=n, sex=sex)
            }
        }
        cat("\n")
    }

    for (class in c("gene", "variant")) {
        for (phenotype in phenotypes)
        {
            cat(paste0(phenotype, ": "))
            subfiles <- grep(paste0(".*chr[0-9X]+_(split_[0-9]+_)*", phenotype, "_(AFR|AMR|EAS|EUR|SAS).*"), files[[class]], value=TRUE)
            pops <- unique(gsub(".*chr[0-9X]+_(split_[0-9]+_)*(.*)_(AFR|AMR|EAS|EUR|SAS).*", "\\3", subfiles))

            for (pop in pops) {
                cat(pop, "\n")
                subsubfiles <- grep(paste0("_", pop, "[_]?[\\.]?"), subfiles, value=TRUE)
                # Ensure that all chromosomes are present
                to_match <- paste0("chr",c(seq(1,22), "X"))
                matched <- to_match %in% unique(gsub(".*(chr[0-9X]+).*", "\\1", subsubfiles))
                if (all(matched))
                {
                    cat(paste0("for population labelling ", pop, " "))
                    cat(paste0("all chromosomes have results files for the phenotype ", phenotype, "\n"))

                    # Combine the results files
                    # First, ensure that no variant files have crept in
                    if (class == "gene") {
                        subsubfiles <- setdiff(subsubfiles, grep("singleAssoc", subsubfiles, value=TRUE))
                    } else {
                        subsubfiles <- intersect(subsubfiles, grep("singleAssoc", subsubfiles, value=TRUE))
                    }
                    cat("files to combine:\n", paste0(subsubfiles, collapse="\n"), "\n")

                    dt <- rbindlist(lapply(subsubfiles, fread))
                    filename <- ifelse(phenotype_info[[phenotype]][[pop]][['binary']],
                        determine_binary_filename(
                            biobank,
                            last_name,
                            analysis_name,
                            phenotype,
                            phenotype_info[[phenotype]][[pop]][['sex']],
                            pop,
                            phenotype_info[[phenotype]][[pop]][['N_case']],
                            phenotype_info[[phenotype]][[pop]][['N_ctrl']],
                            class,
                            date,
                            method,
                            freeze_number),
                        determine_cts_filename(
                            biobank,
                            last_name,
                            analysis_name,
                            phenotype,
                            phenotype_info[[phenotype]][[pop]][['sex']],
                            pop,
                            phenotype_info[[phenotype]][[pop]][['N']],
                            class,
                            date,
                            method,
                            freeze_number)
                        )

                    cat(paste0(data_dir, "/", class, "/sept2023/combined/", filename), "\n")
                    fwrite(dt, quote=FALSE, file=paste0(data_dir, "/", class, "/sept2023/combined/", filename), sep="\t")

                } else {
                    cat(paste0("for the phenotype ", phenotype, " and population labelling ", pop,
                        ", chromosome(s) ", paste(to_match[which(matched)], collapse=", "), " have results\n"))
                    cat(paste0("for the phenotype ", phenotype, " and population labelling ", pop,
                        ", chromosome(s) ", paste(to_match[which(!matched)], collapse=", "), " do not have results\n"))
                }
            }
        }
    }

    # Now, we need to rename the phenotypes according to the correct conventions
    system(paste0("Rscript munge_results_files_Group_names.r",
    	" --folder ", data_dir, "/gene/sept2023/combined",
    	" --out_folder ", out_data_dir, "/gene",
    	" --write")
    )

    system(paste0("Rscript munge_results_files_Group_names.r",
        " --folder ", data_dir, "/variant/sept2023/combined",
        " --out_folder ", out_data_dir, "/variant",
        " --type variant",
        " --write")
    )

    # Finally, upload the cleaned version of the results to the allocated google bucket
} else {
    system(paste0(
        "gsutil -m cp gs://brava-meta-pilot-analysis/pilot-traits-nov-2024-ashg-meta-analysis/",
        biobank, "/gene/* ", out_data_dir, "/gene/")
    )
}

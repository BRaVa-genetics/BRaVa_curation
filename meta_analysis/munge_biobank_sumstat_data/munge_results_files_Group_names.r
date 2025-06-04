#!/bin/Rscript
library(data.table)
library(dplyr)
library(argparse)
library(rlang)

source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")
options(warn = 2)
main <- function(args)
{
	folder_to_check <- args$folder
	type <- args$type
	files_to_check <- grep(paste0("\\.", type, "\\."), dir(folder_to_check, full.names=TRUE), value=TRUE)
	files_to_check <- setdiff(files_to_check, grep("\\.cleaned\\.", files_to_check, value=TRUE))
	cat(paste0("files to check:", paste0(files_to_check, collapse="\n"), "\n"))

	for (f in files_to_check)
	{
		cat(paste0("checking file...", f, "...\n"))
		folder <- gsub("(.*\\/)(.*)", "\\1", f)
		file <- gsub("(.*\\/)(.*)", "\\2", f)
		file_info <- extract_file_info(file)
		cat(paste0("file information:\n", paste0(names(file_info), ": ", file_info, collapse="\n"), "\n\n"))
		renamed <- FALSE
		to_check <- c("dataset", "phenotype", "ancestry", "sex")

		for (check in to_check) {
			if (!(file_info[[check]] %in% names(file_check_information[[check]]))) {
				rename_attempt <- TRUE
	            cat(paste0("attempting to rename ", file_info[[check]], "\n"))
	            for (n in names(file_check_information[[check]])) {
					if (file_info[[check]] %in% file_check_information[[check]][[n]]) {
						cat(paste("Renamed:", file_info[[check]], "->", n, "\n"))
						file_info[[check]] <- n
						renamed <- TRUE
						break
					}
				}
			} else {
				rename_attempt <- FALSE
				cat(paste0("No renaming required for ", file_info[[check]], "\n"))
			}

			if (rename_attempt & !renamed) {
				cat(paste0("Cannot find the appropriate renaming for ", file_info[[check]], "\n"))
			}
		}

        if (file_info$type != type) {
            cat(paste0("This checking script is expecting ", type, " level results\n"))
            stop(paste0(file_info$type, " is not '", type, "'"))
        }

        cat("\nnaming convention checking complete!\n")

        if (!file_info$gz) {
            cat("\nThe file does not appear to be compressed\n",
                "consider compressing before uploading to the cloud\n")
        }

        # Copy the file to a new location with the new name and run all the tests on the
        # proposed filename.
        if (file_info$binary) {
        	required_info <- c("dataset", "last_name", "analysis_name",
        		"phenotype", "freeze_number",  "sex", "ancestry",
        		"n_cases", "n_controls", "software", "type", "date", "split")
        } else {
        	required_info <- c("dataset", "last_name", "analysis_name",
        		"phenotype", "freeze_number",  "sex", "ancestry",
        		"n", "software", "type", "date", "split")
        }

        # Print the old filename
        if (file_info$gz) {
        	file_info$gz <- "gz"
        	required_info <- c(required_info, "gz")
        }

        new_f <- paste0(folder,
        	paste(file_info[required_info], collapse="."))
        cat(paste0(f, " ->\n ", new_f))

        # Now run all the tests on the proposed file names
		exit_status <- system(
			paste0("Rscript ~/Repositories/BRaVa_curation/gene_results_file_checks.r",
				" --just_filename",
				" --file_paths ", new_f,
				" --type ", type)
			)
		
		if (exit_status == 0) {
			cat("successful, rename\n")
			file.rename(f, new_f)
			cat("continue with further checks...\n")

			if (type == "gene") {
				dt <- fread(new_f)
				for (n in names(renaming_header_list)) {
					if (!(n %in% names(dt))) {
						cat(paste("attempting to rename to", n,"\n"))
						if (sum(renaming_header_list[[n]] %in% names(dt)) == 1) {
							to_rename <- which(names(dt) %in% renaming_header_list[[n]])
							cat(paste("renamed:", names(dt)[to_rename], "->", n, "\n"))
							names(dt)[to_rename] <- n
						} else {
							cat(paste("cannot find the column to rename to", n, "\n"))
						}
					}
				}

				for (n in names(renaming_group_list)) {
					to_rename <- which(dt$Group %in% renaming_group_list[[n]])
					if (length(to_rename) > 0) {
						cat("Groups (annotation names) to be renamed...\n")
						cat(paste0(dt$Group[to_rename][1], " -> ", n, "\n"))
						dt$Group[to_rename] <- n
					}
				}

				# Check to ensure that all of the gene names are ensembl IDs
				if (all(grepl("ENSG", dt$Region))) {
					cat("all region names contain ensembl gene IDs\n")
					if (all(grepl("_", dt$Region))) {
						cat("looks like the user has combined ensembl gene ID and genesymbol\n")
						dt[, Region:=gsub(".*_(ENSG[0-9]*)$", "\\1", Region)]
						dt[, Region:=gsub("^(ENSG[0-9]*)_.*", "\\1", Region)]
					}
				} else {
					cat("attempting to resolve to ensembl gene ID\n")
					dt_hgnc <- fread("230117_hgncid_ensembl.txt.gz", select = c("ensembl_gene_id", "hgnc_symbol"))
					dt_hgnc[, Region:=hgnc_symbol]
					dt_hgnc$Region <- ifelse(dt_hgnc$Region == "", dt_hgnc$ensembl_gene_id, dt_hgnc$Region)
					dt_hgnc[, hgnc_symbol:=NULL]
					setkey(dt_hgnc, "Region")
					setkey(dt)
					dt <- merge(dt, dt_hgnc, all.x=TRUE)
					dt[, Region:=ensembl_gene_id]
					dt <- dt %>% filter(Region != "")
				}

				# Remove all other columns
				dt <- dt %>% select(intersect(names(dt), names(renaming_header_list)))
				class <- ifelse(file_info$binary, "binary", "continuous")

				if (any(!(default_gene_result_columns[[class]] %in% names(dt)))) {
					cat("expected columns are missing!\n")
					cat(paste0(setdiff(default_gene_result_columns[[class]], names(dt)), collapse='\n'))
					cat('\n')
					if (any(!(default_gene_result_columns[["minimal"]] %in% names(dt)))) {
						stop("missing a crucial column for meta-analysis")
					} else {
						cat("it's ok though, they're not required for meta-analysis!\n")
					}
				} else {
					cat("all required columns are present!\n")
				}

				if (all((unique(dt$Group) %in% correct_names) | is.na(unique(dt$Group)))) {
					cat("all annotations have the correct names...\n")
					# Include the Cauchy results if the rows are not present
					if (!("Cauchy" %in% unique(dt$Group)))
					{
						cat("Cauchy pvalues are missing, evaluate them...\n")
						dt_cauchy <- dt %>% filter(Group != "Cauchy", !is.na(Pvalue)) %>%
							group_by(Region) %>% summarise(
							Pvalue = cauchy_combination(Pvalue),
							Pvalue_Burden = cauchy_combination(Pvalue_Burden),
							Pvalue_SKAT = cauchy_combination(Pvalue_SKAT),
							number_of_pvals = n()
						) %>% mutate(
							Pvalue = ifelse(Pvalue > 1e+15,
								(1/Pvalue) / pi,
								pcauchy(Pvalue, lower.tail=FALSE)
							),
							Pvalue_Burden = ifelse(Pvalue_Burden > 1e+15,
								(1/Pvalue_Burden) / pi,
								pcauchy(Pvalue_Burden, lower.tail=FALSE)
							),
							Pvalue_SKAT = ifelse(Pvalue_SKAT > 1e+15,
								(1/Pvalue_SKAT) / pi,
								pcauchy(Pvalue_SKAT, lower.tail=FALSE)
							)
						) %>% mutate(
							Pvalue = ifelse(Pvalue > (1 - 1e-10),
								(1 - 1/number_of_pvals), Pvalue
							),
							Pvalue_Burden = ifelse(Pvalue_Burden > (1 - 1e-10),
								(1 - 1/number_of_pvals), Pvalue_Burden
							),
							Pvalue_SKAT = ifelse(Pvalue_SKAT > (1 - 1e-10),
								(1 - 1/number_of_pvals), Pvalue_SKAT
							)
						) %>% select(-number_of_pvals) %>% mutate(Group = "Cauchy")
						dt <- rbind(dt, dt_cauchy, fill=TRUE)
					}

					out_f <- gsub(folder, paste0(args$out_folder, "/"), new_f)
					if (args$write) {
						fwrite(dt, file=out_f, sep="\t", quote=FALSE)
						cat("file written to:\n")
						cat(paste0(out_f, "\n"))
					}
				}
			}

			if (type == "variant")
			{
				cat(paste0("reading in: ", new_f, "\n"))
				dt <- fread(new_f, nrows=1000)
				rename <- FALSE
				marker_fix <- FALSE
				rename_rsid <- FALSE
				rename_chr <- FALSE

				for (n in names(renaming_variant_header_list)) {
					if (!(n %in% names(dt))) {
						cat(paste("attempting to rename to", n,"\n"))
						if (sum(renaming_variant_header_list[[n]] %in% names(dt)) == 1) {
							to_rename <- which(names(dt) %in% renaming_variant_header_list[[n]])
							cat(paste("renamed:", names(dt)[to_rename], "->", n, "\n"))
							names(dt)[to_rename] <- n
							rename <- TRUE
						} else {
							cat(paste("cannot find the column to rename to", n, "\n"))
						}
					}
				}

				# Remove all instances of the ultra-rare variants of various classes in the 
				# genes. We can consider those in separate tests if we need to, and the 
				# remants of them within the constituent variant files can be removed at the 
				# end if need be, so we don't need to worry about them being retained.
				dt <- dt %>% filter(CHR != "UR", !grepl("\\*", MarkerID))

				# Check CHR naming
				if (!all(grepl("^chr[0-9X]+", dt[['CHR']]))) {
					cat("chromosome naming does not match the expected format...\n")
					if (all(grepl("^[0-9XYMT]+$", dt[['CHR']]))) {
						dt[['CHR']] <- paste0('chr', dt[['CHR']])
						rename <- TRUE
						rename_chr <- TRUE
						cat("successfully renamed chromosome column\n")
					} else {
						stop("chromosomes could not be renamed")
					}
				}

				if (!all(grepl("^chr[0-9X]*:[0-9]+:[ACGT]+:[ACGT]+$", dt[['MarkerID']])))
				{
					cat(paste0("marker ID does not match the expected format...\n"))
					cat(paste0("testing to see if it matches other likely formats...\n"))

					# Determine if any of the variants are listed as rsids. If they are, then we convert to
					# chr:pos:ref:alt
					if (any(grepl("^rs[0-9]+", dt[['MarkerID']]))) {
						rename_rsid <- TRUE
					} else if (all(grepl(
						"^chr[0-9X]+[:,/\\_][0-9]+[:,/\\_][ACGT]+[:,/\\_][ACGT]+$",
						dt[['MarkerID']]))) {
						cat(paste0("it matches an alternative potential format, attempting to fix\n"))
						rename <- TRUE
						marker_fix <- TRUE
					} else if (all(grepl("^[0-9X]+[:,/\\_][0-9]+[:,/\\_][ACGT]+[:,/\\_][ACGT]+$",
						dt[['MarkerID']]))) {
						cat(paste0("it matches an alternative potential format, attempting to fix\n"))
						rename <- TRUE
						marker_fix <- TRUE
						dt[['MarkerID']] <- paste0("chr", dt[['MarkerID']])
					} else {
						stop("cannot determine the correct renaming of the variant IDs, please check the file")
					}
				}

				class <- ifelse(file_info$binary, "binary", "continuous")

				if (any(!(default_variant_result_columns[[class]] %in% names(dt)))) {
					cat("expected columns are missing!\n")
					cat(paste0(setdiff(default_variant_result_columns[[class]], names(dt)), collapse='\n'))
					cat('\n')
					if (any(!(default_variant_result_columns[["minimal"]] %in% names(dt)))) {
						stop("missing a crucial column for meta-analysis")
					} else {
						cat("it's ok though, they're not required for meta-analysis!\n")
					}
				} else {
					cat("all required columns are present!\n")
				}

				out_f <- gsub(folder, paste0(args$out_folder, "/"), new_f)
				if (args$write) {
					cat('writing the file...\n')
					if (rename_rsid) {
						cat("running complete rename, using CHR:POS:A1:A2 as the MarkerID\n")
						new_names <- names(dt)
						dt <- fread(new_f)
						cat(paste0("total number of variants is: ", nrow(dt), "\n"))
						names(dt) <- new_names
						dt <- dt %>% filter(CHR != "UR", !grepl("\\*", MarkerID))
						if (rename_chr) {
							dt[['CHR']] <- paste0('chr', dt[['CHR']])
						}
						cat(paste0("following removal of ultra-rares and strangely code variants,",
							" total number of variants is: ", nrow(dt), "\n"))
						dt <- data.table(dt)
						dt[, MarkerID := paste(CHR, POS, Allele1, Allele2, sep=":")]
						# Remove all other columns
						dt <- dt %>% select(intersect(names(dt), names(renaming_variant_header_list)))
						dt <- data.table(dt)
						fwrite(dt, file=out_f, sep="\t", quote=FALSE)
						cat("file written to:\n")
						cat(paste0(out_f, "\n"))
					} else if (rename) {
						new_names <- names(dt)
						dt <- fread(new_f)
						names(dt) <- new_names
						dt <- dt %>% filter(CHR != "UR", !grepl("\\*", MarkerID))
						if (rename_chr) {
							dt[['CHR']] <- paste0('chr', dt[['CHR']])
						}
						dt <- data.table(dt)
						if (marker_fix) {
							dt[, MarkerID := gsub(
								"^(chr)?([0-9X]*)[:,/\\_]([0-9]+)[:,/\\_]([ACGT]+)[:,/\\_]([ACGT]+)$",
								"chr\\2:\\3:\\4:\\5", MarkerID)]
						}
						# Remove all other columns
						dt <- dt %>% select(intersect(names(dt), names(renaming_variant_header_list)))
						dt <- data.table(dt)
						fwrite(dt, file=out_f, sep="\t", quote=FALSE)
						cat("file written to:\n")
						cat(paste0(out_f, "\n"))
					} else if (out_f != new_f) {
						cat("just need to rename the file:\n")
						cat(paste0(out_f, "\n"))
						file.rename(new_f, out_f)
					} else {
						cat("nothing required, the file is ready to go:\n")
						cat(paste0(out_f, "\n"))
					}
				}
			}

		} else {
			warning("Error: file name is not compatible with the requested format")
		}
	}
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--folder", default=NULL, required=TRUE,
    help="folder to check")
parser$add_argument("--out_folder", default=NULL, required=TRUE,
    help="folder to export results to")
parser$add_argument("--write", default=FALSE, action="store_true",
    help="write the resultant file?")
parser$add_argument("--type", default="gene", required=FALSE)
args <- parser$parse_args()

main(args)


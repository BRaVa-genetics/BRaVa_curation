#!/bin/Rscript
library(data.table)
library(dplyr)
library(stringr)
library(argparse)

source("../meta_analysis_utils.r")

# For naming files
biobank <- "uk-biobank"
last_name <- "palmer"
analysis_name <- "pilot"
freeze_number <- "JULY23Freeze"
date <- "20240110"
method <- "SAIGE"

dx_data_dir <- "brava/outputs/step2/sept2023/"

data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

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


# For naming files
biobank <- "uk-biobank"
last_name <- "baya"
analysis_name <- "pilot"
freeze_number <- "JULY23Freeze"
date <- "20250606"
method <- "regenie"
pop <- "EUR"
phenotype <- "Height"
sex <- "ALL"

# Next, deal with height for EUR - note that this did not converge in SAIGE
dx_data_var <- "/nbaya/regenie/data/step2/variant_tests/Height/EUR/regenie_variant_test.EUR.Height.chr*.regenie.gz"
dx_data_gene <- "/nbaya/regenie/data/step2/group_tests/Height/EUR/*.gz"
dx_data_annot <- "/nbaya/regenie/data/annotations/v7/*annotations*"
data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/raw/height")
out_data_dir <- paste0("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks/", biobank, "/cleaned")

system(paste0("mkdir -p ", data_dir, "/gene/june2025/combined"))
system(paste0("mkdir -p ", data_dir, "/variant/june2025/combined"))
system(paste0("mkdir -p ", data_dir, "/annotations/june2025/combined"))

if (download) {
    # Replace with dx download, and we now need to combine across chromosomes
    # Also, we should flag if a chromosome is missing!
    system(paste0("dx download ", dx_data_var, " -o ", data_dir , "/variant/june2025/"))
    system(paste0("dx download ", dx_data_gene, " -o ", data_dir , "/gene/june2025/"))
    system(paste0("dx download ", dx_data_annot, " -o ", data_dir, "/annotations/june2025/"))
}

# Now munge
# Regenie group test results
class <- "gene"
dt_list <- list()
for (file in dir(paste0(data_dir , "/gene/june2025/"),
    pattern=paste0(".*group.*", pop, ".*gz$"), full.names=TRUE)) {
    dt_list[[file]] <- fread(file, header=TRUE)
}
dt_regenie <- rbindlist(dt_list, fill=TRUE)
N <- dt_regenie$N[1]

dt_regenie <- dt_regenie %>% filter(
    TEST %in% c("ADD", "ADD-SKAT", "ADD-SKATO"),
    !grepl("singleton", ALLELE1))
dt_regenie <- dt_regenie %>% mutate(
    Region = gsub("(ENSG[0-9]+)\\..*", "\\1", ID),
    max_MAF = as.numeric(gsub(".*(0\\.0.*$)", "\\1", ID)),
    Group = gsub("^([^\\.]*)\\..*", "\\1", ALLELE1)
    ) %>% mutate(Group = gsub(":", ";", Group))
dt_regenie <- dt_regenie %>% rename(BETA_Burden=BETA, SE_Burden=SE) %>% 
    select(-c("ALLELE1", "ALLELE0", "ID", "CHROM", "GENPOS", "EXTRA"))

dt_regenie_burden <- data.table(dt_regenie %>% filter(TEST == "ADD") %>% 
    mutate(Pvalue_Burden = 10^(-LOG10P)) %>% select(-c("N", "TEST")) %>% 
    rename(CHISQ_burden=CHISQ, LOG10P_burden=LOG10P))
dt_regenie_SKAT <- data.table(dt_regenie %>% filter(TEST == "ADD-SKAT") %>%
    mutate(Pvalue_SKAT = 10^(-LOG10P)) %>%
    select(-c("A1FREQ", "N", "TEST", "BETA_Burden", "SE_Burden")) %>%
    rename(CHISQ_SKAT=CHISQ, LOG10P_SKAT=LOG10P))
dt_regenie_SKATO <- data.table(dt_regenie %>% filter(TEST == "ADD-SKATO") %>%
    mutate(Pvalue = 10^(-LOG10P)) %>%
    select(-c("A1FREQ", "N", "TEST", "BETA_Burden", "SE_Burden")))

# Get the naming convention right
setkeyv(dt_regenie_burden, c("Region", "Group", "max_MAF"))
setkeyv(dt_regenie_SKAT, c("Region", "Group", "max_MAF"))
setkeyv(dt_regenie_SKATO, c("Region", "Group", "max_MAF"))

dt_regenie <- merge(merge(dt_regenie_burden, dt_regenie_SKAT), dt_regenie_SKATO)
dt_regenie <- data.table(dt_regenie)
setkeyv(dt_regenie, c("Region", "Group", "max_MAF"))
dt_regenie <- dt_regenie %>% select(Region, Group, max_MAF, Pvalue, Pvalue_Burden, Pvalue_SKAT,
        BETA_Burden, SE_Burden)
filename <- determine_cts_filename(
    biobank,
    last_name,
    analysis_name,
    phenotype,
    sex,
    pop,
    N,
    class,
    date,
    method,
    freeze_number)
cat(paste0(data_dir, "/", class, "/june2025/combined/", filename), "\n")
fwrite(dt_regenie, quote=FALSE, file=paste0(data_dir, "/", class, "/june2025/combined/", filename), sep="\t")

# Now the same for the variant based tests
class <- "variant"
dt_list <- list()
for (file in dir(paste0(data_dir , "/variant/june2025/"),
    pattern=paste0(pop, ".*gz$"), full.names=TRUE)) {
    dt_list[[file]] <- fread(file, header=TRUE)
}
dt_regenie <- rbindlist(dt_list, fill=TRUE)
dt_regenie <- dt_regenie %>% rename(
    CHR=CHROM, POS=GENPOS, MarkerID=ID,
    Allele1=ALLELE0, Allele2=ALLELE1, AF_Allele2=A1FREQ)
N_tot <- N
dt_regenie <- dt_regenie %>% 
    mutate(
        AC_Allele2=round(2*N*AF_Allele2),
        MissingRate=N/N_tot) %>% 
    mutate(
        AF_Allele2=AC_Allele2/(2*N_tot),
        N=N_tot,
        Tstat=CHISQ/BETA) %>%
    mutate(var=Tstat^2/CHISQ,
        `p.value` = 10^(-LOG10P)) %>%
    select(CHR, POS, MarkerID, Allele1, Allele2,
        AC_Allele2, AF_Allele2, MissingRate, BETA,
        SE, Tstat, var, p.value, N)

filename <- determine_cts_filename(
    biobank,
    last_name,
    analysis_name,
    phenotype,
    sex,
    pop,
    N,
    class,
    date,
    method,
    freeze_number)
fwrite(dt_regenie, quote=FALSE, file=paste0(data_dir, "/", class, "/june2025/combined/", filename), sep="\t")
dt_regenie <- fread(paste0(data_dir, "/", class, "/june2025/combined/", filename))
# Also need the annotation files as well
# Read in the annotation file and merge
dt_list <- list()
for (file in dir(paste0(data_dir, "/annotations/june2025"),
        pattern=".txt", full.names=TRUE)) {
    dt_list[[file]] <- fread(file, header=FALSE)
}
dt_annot <- rbindlist(dt_list, fill=TRUE) %>% rename(MarkerID=V1, Region=V2, Group=V3)
dt_annot <- data.table(dt_annot)

dt_regenie <- data.table(dt_regenie)
setkeyv(dt_regenie, "MarkerID")
setkeyv(dt_annot, "MarkerID")
dt_regenie <- merge(dt_regenie, dt_annot) %>% filter(Group != "non_coding")

compute_burden <- function(thresh) {
    dt_gene <- dt_regenie %>% filter(
            pmin(AF_Allele2, 1-AF_Allele2) < thresh,
            Group != "non_coding") %>% 
        group_by(Region, Group) %>% mutate(w=dbeta(pmin(AF_Allele2, 1-AF_Allele2), 1, 25)) %>%
        summarise(
            BETA_Burden = sum(Tstat*w)/sum(var*w^2),
            SE_Burden=sqrt(1/sum(var*w^2)),
            .groups = "drop")

    # Finally determine the union of damaging missense and pLoF, and of the lot
    dt_gene_dm_pLoF <- dt_regenie %>% filter(pmin(AF_Allele2, 1-AF_Allele2) < thresh,
            Group %in% c("damaging_missense_or_protein_altering", "pLoF")) %>% 
        group_by(Region) %>% mutate(w=dbeta(pmin(AF_Allele2, 1-AF_Allele2), 1, 25)) %>%
        summarise(
            BETA_Burden = sum(Tstat*w)/sum(var*w^2),
            SE_Burden=sqrt(1/sum(var*w^2)),
            .groups = "drop") %>% 
        mutate(Group = "pLoF;damaging_missense_or_protein_altering")

    dt_gene_any <- dt_regenie %>% filter(pmin(AF_Allele2, 1-AF_Allele2) < thresh,
            Group %in% c(
                "pLoF",
                "damaging_missense_or_protein_altering",
                "other_missense_or_protein_altering",
                "synonymous")) %>% 
        group_by(Region) %>% mutate(w=dbeta(pmin(AF_Allele2, 1-AF_Allele2), 1, 25)) %>%
        summarise(
            BETA_Burden = sum(Tstat*w)/sum(var*w^2),
            SE_Burden=sqrt(1/sum(var*w^2)),
            .groups = "drop") %>% 
        mutate(Group = "pLoF;damaging_missense_or_protein_altering;other_missense_or_protein_altering;synonymous")
    return(rbind(dt_gene, dt_gene_dm_pLoF, dt_gene_any) %>% mutate(max_MAF=thresh))
}

dt_list <- list()
for (t in c(0.01, 0.001, 0.0001)) {
    cat(t, "\n")
    dt_list[[as.character(t)]] <- compute_burden(t)
}
dt_gene <- rbindlist(dt_list)
class <- "gene"
filename <- determine_cts_filename(
    biobank,
    last_name,
    analysis_name,
    phenotype,
    sex,
    pop,
    N,
    class,
    date,
    method,
    freeze_number)
fwrite(dt_gene, paste0(data_dir, "/", class, "/june2025/combined/munged.", filename), sep='\t', quote=FALSE)

dt_gene <- fread(paste0(data_dir, "/", class, "/june2025/combined/munged.", filename))
dt_gene_regenie <- fread(paste0(data_dir, "/", class, "/june2025/combined/", filename))

setkeyv(dt_gene, c("Region", "Group", "max_MAF"))
setkeyv(dt_gene_regenie, c("Region", "Group", "max_MAF"))
dt_gene <- merge(dt_gene, dt_gene_regenie %>% select(-c("BETA_Burden", "SE_Burden")))
fwrite(dt_gene, paste0(data_dir, "/", class, "/june2025/combined/", filename))
system(paste("rm", paste0(data_dir, "/", class, "/june2025/combined/munged.", filename)))

# Now, we need to rename the phenotypes according to the correct conventions
system(paste0("Rscript munge_results_files_Group_names.r",
    " --folder ", data_dir, "/gene/june2025/combined",
    " --out_folder ", out_data_dir, "/gene",
    " --write")
)

system(paste0("Rscript munge_results_files_Group_names.r",
    " --folder ", data_dir, "/variant/june2025/combined",
    " --out_folder ", out_data_dir, "/variant",
    " --type variant",
    " --write")
)

# pdf(width=10, height=3, file="maf.pdf")
# p <- ggplot(dt_gene, aes(x=BETA_Burden.x, y=BETA_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=SE_Burden.x, y=SE_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=(BETA_Burden.x/SE_Burden.x), y=(BETA_Burden.y/SE_Burden.y))) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# dev.off()

# pdf(width=20, height=10, file="group.pdf")
# p <- ggplot(dt_gene, aes(x=BETA_Burden.x, y=BETA_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=SE_Burden.x, y=SE_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=(BETA_Burden.x/SE_Burden.x), y=(BETA_Burden.y/SE_Burden.y))) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# dev.off()

# dt_gene <- dt_gene %>% filter(Pvalue_Burden < 0.01)
# pdf(width=10, height=3, file="maf_p.pdf")
# p <- ggplot(dt_gene, aes(x=BETA_Burden.x, y=BETA_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=SE_Burden.x, y=SE_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=(BETA_Burden.x/SE_Burden.x), y=(BETA_Burden.y/SE_Burden.y))) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~max_MAF) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# dev.off()

# pdf(width=20, height=10, file="group_p.pdf")
# p <- ggplot(dt_gene, aes(x=BETA_Burden.x, y=BETA_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=SE_Burden.x, y=SE_Burden.y)) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# p <- ggplot(dt_gene, aes(x=(BETA_Burden.x/SE_Burden.x), y=(BETA_Burden.y/SE_Burden.y))) + geom_hex(bins = 200) +
#   scale_fill_viridis_c(option = "plasma") + theme_minimal() + facet_wrap(~Group) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")
# print(p)
# dev.off()

#!/usr/bin/env Rscript
library(data.table)
library(biomaRt)
# required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(future.apply)
library(progressr)

source("meta_analysis_utils.r")

# Input directory (same as before)
files <- dir(path="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100", full.names=TRUE)

ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Retrieve Ensembl IDs, start and end positions and external name
gene_info <- getBM(attributes = c(
  "ensembl_gene_id",
  "external_gene_name",
  "chromosome_name",
  "start_position",
  "end_position"),
  mart = ensembl_dataset)

gene_info <- data.table(gene_info)
# Build a named vector for lookup (safer)
gene_map <- setNames(gene_info$external_gene_name, gene_info$ensembl_gene_id)

extract_file_info_meta <- function(filename)
{
  gz <- ifelse(grepl(".gz$", filename), TRUE, FALSE)
  filename_nogz <- gsub(".gz$", "", filename)
  filename_nogz <- gsub("cleaned.", "", filename_nogz)
  file_info <- as.list(strsplit(filename_nogz, split="\\.")[[1]])
  file_info[[1]] <- gsub(".*/", "", file_info[[1]])
  file_info$phenotype <- gsub("^([A-Za-z0-9]+).*", "\\1", file_info[[1]])
  file_info$meta <- ifelse(length(file_info) >= 2 && file_info[[2]] == "tsv", "ALL", ifelse(length(file_info) >= 2, file_info[[2]], NA))
  file_info$gz <- gz
  file_info <- file_info[c("phenotype", "meta", "gz")]
  return(file_info)
}

plan(multisession, workers = 32)   # adjust workers if needed

# --- progress setup ---
handlers("progress")  # nice progress bar

# threshold
P_burden <- 0.05/(20000*3*2)

allowed_groups <- c(
  "pLoF",
  "damaging_missense_or_protein_altering",
  "pLoF;damaging_missense_or_protein_altering"
)

allowed_mafs <- c(0.001, 0.0001)

with_progress({

  p <- progressor(along = files)

  dt_list <- future_lapply(files, function(file) {

    p(message = basename(file))   # <-- progress update

    raw <- tryCatch(fread(file), error = function(e) return(NULL))
    if (is.null(raw)) return(NULL)

    needed_cols <- c(
      "Region", "Pvalue", "Group", "max_MAF",
      "type", "class", "BETA_Burden", "SE_Burden"
    )

    if (!all(needed_cols %in% names(raw))) return(NULL)

    raw[, max_MAF := suppressWarnings(as.numeric(max_MAF))]

    dt <- raw[
      Group %in% allowed_groups &
      max_MAF %in% allowed_mafs &
      (type == "Inverse variance weighted" & class == "Burden" & Pvalue < P_burden),
      .(Region, Group, max_MAF, class, type, Pvalue, BETA_Burden, SE_Burden)
    ]

    if (nrow(dt) == 0) return(NULL)

    file_info <- extract_file_info_meta(file)

    dt[, `:=`(
      phenotype = file_info$phenotype,
      meta = file_info$meta,
      file = file
    )]

    dt
  })
})

dt_final <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Map Region -> external gene name; safe if some Regions not present will become NA
dt_final$external_gene_name <- gene_map[as.character(dt_final$Region)]

# Create named vectors for fast lookup
renaming_map <- unlist(renaming_phenotype_list)
broad_map <- unlist(phenotype_broad_categories)

class_map <- c(
  setNames(rep("binary", length(phenotype_class$binary)), phenotype_class$binary),
  setNames(rep("continuous", length(phenotype_class$continuous)), phenotype_class$continuous)
)

# Add columns to dt_final (NA where mapping missing)
dt_final[, phenotype_full := renaming_map[phenotype]]
dt_final[, phenotype_broad_category := broad_map[phenotype]]
dt_final[, phenotype_class := class_map[phenotype]]

# Setting P_burden and P_skat as earlier gives old results.
fwrite(dt_final, sep='\t', file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_burden_only.tsv.gz")

# --- parallel setup ---
plan(multisession, workers = 32)   # adjust workers if needed

# --- progress setup ---
handlers("progress")  # nice progress bar

# thresholds
P_burden <- 0.05/(20000*3*3*2)
P_skat   <- 0.05/(20000*3*3*2)

allowed_groups <- c(
  "pLoF",
  "damaging_missense_or_protein_altering",
  "pLoF;damaging_missense_or_protein_altering"
)

allowed_mafs <- c(0.001, 0.0001)

with_progress({

  p <- progressor(along = files)

  dt_list <- future_lapply(files, function(file) {

    p(message = basename(file))   # <-- progress update

    raw <- tryCatch(fread(file), error = function(e) return(NULL))
    if (is.null(raw)) return(NULL)

    needed_cols <- c(
      "Region", "Pvalue", "Group", "max_MAF",
      "type", "class", "BETA_Burden", "SE_Burden"
    )

    if (!all(needed_cols %in% names(raw))) return(NULL)

    raw[, max_MAF := suppressWarnings(as.numeric(max_MAF))]

    dt <- raw[
      Group %in% allowed_groups &
      max_MAF %in% allowed_mafs &
      (
        (type == "Inverse variance weighted" & class == "Burden" & Pvalue < P_burden) |
        (type == "Stouffer" & class %in% c("SKAT", "SKAT-O") & Pvalue < P_skat)
      ),
      .(Region, Group, max_MAF, class, type, Pvalue, BETA_Burden, SE_Burden)
    ]

    if (nrow(dt) == 0) return(NULL)

    file_info <- extract_file_info_meta(file)

    dt[, `:=`(
      phenotype = file_info$phenotype,
      meta = file_info$meta,
      file = file
    )]

    dt
  })
})

dt_final <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Map Region -> external gene name; safe if some Regions not present will become NA
dt_final$external_gene_name <- gene_map[as.character(dt_final$Region)]

# Create named vectors for fast lookup
renaming_map <- unlist(renaming_phenotype_list)
broad_map <- unlist(phenotype_broad_categories)

class_map <- c(
  setNames(rep("binary", length(phenotype_class$binary)), phenotype_class$binary),
  setNames(rep("continuous", length(phenotype_class$continuous)), phenotype_class$continuous)
)

# Add columns to dt_final (NA where mapping missing)
dt_final[, phenotype_full := renaming_map[phenotype]]
dt_final[, phenotype_broad_category := broad_map[phenotype]]
dt_final[, phenotype_class := class_map[phenotype]]

# Setting P_burden and P_skat as earlier gives old results.
fwrite(dt_final, sep='\t', file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent.tsv.gz")

with_progress({

  p <- progressor(along = files)

  dt_list <- future_lapply(files, function(file) {

    p(message = basename(file))

    raw <- tryCatch(fread(file), error = function(e) return(NULL))
    if (is.null(raw)) return(NULL)

    needed_cols <- c(
      "Region", "Pvalue", "Group", "max_MAF",
      "type", "class"
    )

    if (!all(needed_cols %in% names(raw))) return(NULL)

    raw[, max_MAF := suppressWarnings(as.numeric(max_MAF))]

    dt <- raw[
      Group %in% allowed_groups &
      max_MAF %in% allowed_mafs &
      (
        (type == "Inverse variance weighted" & class == "Burden") |
        (type == "Stouffer" & class %in% c("SKAT", "SKAT-O"))
      ),
      .(Region, Group, max_MAF, class, type, Pvalue)
    ]
    dt <- dt[!is.na(Pvalue) & is.finite(Pvalue)]

    if (nrow(dt) == 0) return(NULL)

    # --- metadata ---
    file_info <- extract_file_info_meta(file)

    dt[, `:=`(
      phenotype = file_info$phenotype,
      meta = file_info$meta
    )]

    dt[, Pvalue := pmin(pmax(Pvalue, 1e-300), 0.99)]

    # --- Cauchy combination ---
    dt_cauchy <- dt[, {

      stat <- cauchy_combination(Pvalue, weights = rep(1, .N))
      n_p <- .N

      p_meta <- if (stat > 1e+15) {
        (1 / stat) / pi
      } else {
        pcauchy(stat, lower.tail = FALSE)
      }

      if (p_meta > (1 - 1e-10)) {
        p_meta <- (1 - 1 / n_p)
      }

      .(
        Cauchy_stat = stat,
        number_of_pvals = n_p,
        Cauchy_Pvalue = p_meta
      )

    }, by = .(Region, phenotype, meta)]

    dt_cauchy
  })
})

# --- combine all files ---
dt_cauchy_all <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
sig_pairs <- dt_cauchy_all[Cauchy_Pvalue < 0.05/20000]
sig_pairs <- unique(sig_pairs[, .(Region, phenotype, meta)])
sig_split <- split(
  sig_pairs,
  by = c("phenotype", "meta"),
  keep.by = FALSE
)

with_progress({

  p <- progressor(along = files)

  dt_list2 <- future_lapply(files, function(file) {

    p(message = basename(file))

    file_info <- extract_file_info_meta(file)
    pheno <- file_info$phenotype
    meta  <- file_info$meta

    key <- paste(pheno, meta, sep = ".")

    # Skip if no significant pairs for this combo
    if (!key %in% names(sig_split)) return(NULL)

    targets <- sig_split[[key]]$Region

    raw <- tryCatch(fread(file), error = function(e) return(NULL))
    if (is.null(raw)) return(NULL)

    raw[, max_MAF := suppressWarnings(as.numeric(max_MAF))]

    dt <- raw[
      Region %in% targets &
      Group %in% allowed_groups &
      max_MAF %in% allowed_mafs &
      (
        (type == "Inverse variance weighted" & class == "Burden") |
        (type == "Stouffer" & class %in% c("SKAT", "SKAT-O"))
      )
    ]

    if (nrow(dt) == 0) return(NULL)

    dt[, `:=`(
      phenotype = pheno,
      meta = meta,
      file = file
    )]

    dt
  })
})

dt_final_filtered <- rbindlist(dt_list2, use.names = TRUE, fill = TRUE)

# Map Region -> external gene name; safe if some Regions not present will become NA
dt_final_filtered$external_gene_name <- gene_map[as.character(dt_final_filtered$Region)]

# Add columns to dt_final (NA where mapping missing)
dt_final_filtered[, phenotype_full := renaming_map[phenotype]]
dt_final_filtered[, phenotype_broad_category := broad_map[phenotype]]
dt_final_filtered[, phenotype_class := class_map[phenotype]]

# Make sure Region and phenotype ordering is preserved; create unique rows for each (they should already be per-file)
dt_final_filtered <- unique(dt_final_filtered)

# This is the updated results - Cauchy sig genes pulled out.
fwrite(dt_final_filtered, sep='\t', file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_cauchy.tsv.gz")

p_thres <- 0.05/(20000*3*3*2)

summary_dt <- dt_final_filtered[
  , {
      pvals <- as.numeric(Pvalue)
      sig <- !is.na(pvals) & is.finite(pvals) & pvals < p_thres

      min_idx <- which.min(pvals)

      .(
        external_gene_name = unique(na.omit(external_gene_name))[1L],
        phenotype_full = unique(na.omit(phenotype_full))[1L],
        phenotype_broad_category = unique(na.omit(phenotype_broad_category))[1L],
        phenotype_class = unique(na.omit(phenotype_class))[1L],

        metas = paste(sort(unique(meta[sig])), collapse = "; "),
        n_metas = length(unique(meta[sig])),

        groups = paste(sort(unique(Group[sig])), collapse = "; "),
        max_MAFs = paste(sort(unique(as.character(max_MAF[sig]))), collapse = "; "),
        classes = paste(sort(unique(class[sig])), collapse = "; "),
        types = paste(sort(unique(type[sig])), collapse = "; "),

        min_pvalue = min(pvals, na.rm = TRUE),
        meta_min_p = meta[min_idx][1L],
        BETA_Burden_min_p = BETA_Burden[min_idx][1L],
        SE_Burden_min_p = SE_Burden[min_idx][1L],

        best_test = paste0(
          class[min_idx][1L], ":", type[min_idx][1L], ":",
          Group[min_idx][1L], ":maf=", max_MAF[min_idx][1L]
        ),

        any_burden_signif = any(class == "Burden" & sig),
        any_skat_signif = any(class %in% c("SKAT", "SKAT-O") & sig),

        metas_burden = paste(sort(unique(meta[class == "Burden" & sig])), collapse = "; "),
        metas_skat = paste(sort(unique(meta[class %in% c("SKAT", "SKAT-O") & sig])), collapse = "; "),

        n_hits = sum(sig) 
      )
    },
  by = .(Region, phenotype)
]

# Now, determine the most signficant Cauchy combination Pvalue (where at least one is significant)

summary_dt_cauchy <- dt_cauchy_all[Cauchy_Pvalue < 0.05/20000][
  , .(
    metas_cauchy = paste(sort(unique(meta)), collapse = "; "),
    n_metas_cauchy = length(unique(meta)),
    min_cauchy_pvalue = min(as.numeric(Cauchy_Pvalue), na.rm = TRUE),
    best_cauchy = meta[which.min(as.numeric(Cauchy_Pvalue))][1L],
    n_cauchy_hits = .N
    ),
  by = .(Region, phenotype)
]

# Map Region -> external gene name; safe if some Regions not present will become NA
summary_dt_cauchy$external_gene_name <- gene_map[as.character(summary_dt_cauchy$Region)]

# Add columns to dt_final (NA where mapping missing)
summary_dt_cauchy[, phenotype_full := renaming_map[phenotype]]
summary_dt_cauchy[, phenotype_broad_category := broad_map[phenotype]]
summary_dt_cauchy[, phenotype_class := class_map[phenotype]]

# Finally, merge in the Cauchy Pvalue
summary_dt <- summary_dt[
  summary_dt_cauchy,
  on = .(Region, phenotype, external_gene_name, phenotype_full, phenotype_broad_category, phenotype_class)
]

# If external_gene_name missing, fill with Region
summary_dt[is.na(external_gene_name) | external_gene_name == "", external_gene_name := Region]

# Optionally reorder columns to something readable
setcolorder(summary_dt, c("Region", "external_gene_name", "phenotype", "phenotype_full",
                          "phenotype_broad_category", "phenotype_class",
                          "n_hits", "n_metas", "metas", "metas_burden", "metas_skat",
                          "groups", "max_MAFs", "classes", "types",
                          "min_pvalue", "meta_min_p", "BETA_Burden_min_p", "SE_Burden_min_p",
                          "best_test", "any_burden_signif", "any_skat_signif",
                          "n_cauchy_hits", "n_metas_cauchy", "metas_cauchy",
                          "min_cauchy_pvalue", "best_cauchy"))

fwrite(summary_dt, sep='\t', file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz")

summary_dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz")
summary_dt_for_agent <- unique(summary_dt %>% select(Region, external_gene_name, phenotype, phenotype_full, min_cauchy_pvalue))
# This is the file that should be passed to the AI agent
fwrite(summary_dt_for_agent, sep="\t", file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_final_summary.tsv.gz")

# Preprocess once - this is prepping to include all the results for biobank only (gene, phenotypes)
# if they are significant in meta - and get the results for all the other meta-analyses.
dt_final <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_cauchy.tsv.gz")
dt_final <- unique(dt_final[, .(phenotype, Region, Group, max_MAF)])

# Split once
dt_final_split <- split(dt_final, by = "phenotype", keep.by = FALSE)

# File list
files <- dir(
  path = "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene/n_cases_100",
  full.names = TRUE,
  recursive = TRUE
)
files <- files[!grepl("/minus_|/just_uk-biobank_and_all-of-us", files)]

# Progress-aware parallel loop
with_progress({
  p <- progressor(along = files)

  dt_list <- future_lapply(files, function(file) {

    p(file)  # update progress

    file_info <- extract_file_info_meta(file)
    pheno <- file_info$phenotype

    if (!pheno %in% names(dt_final_split)) return(NULL)

    dt_tmp <- dt_final_split[[pheno]]
    setkey(dt_tmp, Region, Group, max_MAF)
    dt_tmp <- unique(dt_tmp, by = key(dt_tmp))

    raw <- fread(file)
    raw[, max_MAF := as.numeric(max_MAF)]
    setkey(raw, Region, Group, max_MAF)

    dt <- raw[dt_tmp, nomatch = 0]

    dt <- dt[
      (type == "Inverse variance weighted" & class == "Burden") |
      (type == "Stouffer" & class %in% c("SKAT", "SKAT-O"))
    ]

    if (nrow(dt) == 0) return(NULL)

    dt <- dt[, .(
      Region, Group, max_MAF, class, type,
      Pvalue, BETA_Burden, SE_Burden
    )]

    dt[, `:=`(
      phenotype = pheno,
      meta = file_info$meta,
      file = file
    )]

    dt
  }, future.seed = TRUE)
})

# Combine
dt_final_superset <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

# Map Region -> external gene name; safe if some Regions not present will become NA
dt_final_superset$external_gene_name <- gene_map[as.character(dt_final_superset$Region)]

# Add columns to dt_final (NA where mapping missing)
dt_final_superset[, phenotype_full := renaming_map[phenotype]]
dt_final_superset[, phenotype_broad_category := broad_map[phenotype]]
dt_final_superset[, phenotype_class := class_map[phenotype]]

# Make sure Region and phenotype ordering is preserved; create unique rows for each (they should already be per-file)
dt_final_superset <- unique(dt_final_superset)
fwrite(dt_final_superset, sep='\t', file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_superset.tsv.gz")
dt_final_superset <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_superset.tsv.gz")

# Read in the results sent back by Jeremy
# dt <- fread("manuscript_figures/novelty_agent_results.tsv")
# Uploaded agent results to the cluster
dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/novelty_agent_results.tsv")
summary_dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz") %>%
  rename(gene=external_gene_name, phenotype_ID=phenotype) %>% rename(phenotype=phenotype_full)

# Determine whether additional genes need to be included.
# If so, run and reupload the appended version.
# This is now done (hence commented below)
# dt_new <- setdiff(unique(summary_dt %>% select(gene, phenotype)), dt %>% select(gene, phenotype))
# setkeyv(dt_new, c("gene", "phenotype"))
# setkeyv(summary_dt, c("gene", "phenotype"))
# fwrite(file="/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/missed_for_agent.tsv", merge(dt_new, summary_dt) %>% select(gene, phenotype, min_pvalue))
# # To merge any new things in:
# dt <- fread("novelty_agent_results.tsv")
# dt_new <- fread("new_novelty_agent_results.tsv") %>% 
#   rename(literature_verdict = verdict_literature,
#     open_targets_verdict = verdict_open_targets,
#     max_verdict = verdict_max)
# dt <- rbind(dt, dt_new, use.names=TRUE)
# fwrite(dt, file="novelty_agent_results.tsv", sep='\t', quote=FALSE)

# Want further information detailing the most significant association in a given (ancestry biobank) pair
# Keen to determine which associations would not have been found for a single biobank.
setkeyv(dt, c("gene", "phenotype"))
setkeyv(summary_dt, c("gene", "phenotype"))

dt <- merge(dt, summary_dt)
dt <- dt %>% select(gene, Region, phenotype, phenotype_ID,  phenotype_broad_category,
  literature_verdict, open_targets_verdict, max_verdict, justification, citations,
  phenotype_class, n_hits, n_metas, metas_burden, metas_skat, metas, groups, max_MAFs,
  classes, types, min_pvalue, meta_min_p, best_test, any_burden_signif, any_skat_signif,
  n_cauchy_hits, metas_cauchy, min_cauchy_pvalue, best_cauchy)

dt <- dt %>% rename(
  `Gene symbol` = gene,
  `Gene ID` = Region,
  Description = phenotype,
  `Phenotype ID` = phenotype_ID,
  `Phenotype broad category` = phenotype_broad_category,
  `Literature verdict` = literature_verdict,
  `Open targets verdict` = open_targets_verdict,
  `Max verdict` = max_verdict,
  Justification = justification,
  `Citations (PubMed IDs)` = citations,
  `Binary or continuous` = phenotype_class,
  `N significant associations across (max MAF, mask, meta-analyses)` = n_hits,
  `N significant meta-analyses` = n_metas,
  `Significant meta-analyses` = metas_burden,
  `Significant burden meta-analyses` = metas,
  `Significant SKAT/SKAT-O meta-analyses` = metas_skat,
  `Significant masks` = groups,
  `Significant max MAFs` = max_MAFs,
  `Significant classes` = classes,
  `Meta-analysis method` = types,
  `Min P-value of significant associations` = min_pvalue,
  `Meta-analysis with minimum P-value` = meta_min_p,
  `Mask for most significant association` = best_test,
  `Any burden significant` = any_burden_signif,
  `Any SKAT/SKAT-O significant` = any_skat_signif,
  `N significant Cauchy associations across meta-analyses` = n_cauchy_hits,
  `Significant Cauchy` = metas_cauchy,
  `Min P-value of significant Cauchy associations` = min_cauchy_pvalue,
  `Cauchy combination with minimum P-value` = best_cauchy
  )

# Merged version.
fwrite(dt, file = "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz", sep='\t', quote=FALSE)

# Write the code to determine association P-values for all of the significant associations in the gene_phenotype pairs above

# Here, we must remove any files that have been deemed to be inflated.
inflation_file <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz"
dt_inflation <- fread(inflation_file)
dt_inflation <- unique(dt_inflation %>% filter(Group == "synonymous") %>% 
  filter(max_MAF != 0.01,
    lambda_value > 1.3,
    !(lambda_type %in% c("lambda_50_Burden", "lambda_50_SKAT", "lambda_50"))) %>% 
  dplyr::select(phenotype, dataset, ancestry, sex))
# Manual curation, adding the following (biobank, trait) tuples containing spurious 
# associations
dt_inflation <- rbind(dt_inflation, data.table(
  phenotype = c("ColonRectCanc", "Height"),
  dataset = c("egcut", "mgbb"),
  ancestry = c("EUR", "AMR"),
  sex = c("ALL", "ALL"))
) %>% rename(phenotypeID = phenotype, biobank = dataset, pop = ancestry)

# Main lookup table
dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz")
setDT(dt)

# phenotype -> Region lookup, built once
dt_map <- unique(dt[, .(
  phenotype = `Phenotype ID`,
  Region = `Gene ID`
)])
phenotype_regions <- split(dt_map$Region, dt_map$phenotype)

# Inflated-file lookup, built once
inflated_dt <- unique(dt_inflation[, .(
  phenotypeID,
  biobank,
  pop
)])
setkey(inflated_dt, phenotypeID, biobank, pop)

keep_groups <- c(
  "pLoF",
  "damaging_missense_or_protein_altering",
  "pLoF;damaging_missense_or_protein_altering"
)

# Build a flat file table: one row per file, with its biobank
biobanks <- dir("/well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/biobanks", full.names = TRUE)

file_dt <- rbindlist(lapply(biobanks, function(biob) {
  files <- dir(file.path(biob, "cleaned", "gene"), full.names = TRUE)
  if (length(files) == 0L) return(NULL)
  data.table(
    file = files,
    biobank = basename(biob)
  )
}), use.names = TRUE, fill = TRUE)

# Parallel loop with progress
with_progress({
  p <- progressor(along = file_dt$file)

  dt_list <- future_lapply(seq_len(nrow(file_dt)), function(i) {
    p(basename(file_dt$file[i]))

    file <- file_dt$file[i]
    biobank <- file_dt$biobank[i]

    info <- extract_file_info(file)
    pheno <- info$phenotype
    ancestry <- info$ancestry

    # Skip inflated sumstats
    if (inflated_dt[.(pheno, biobank, ancestry), nomatch = 0L, .N] > 0L) {
      return(NULL)
    }

    regions <- phenotype_regions[[pheno]]
    if (is.null(regions) || length(regions) == 0L) return(NULL)

    dt_tmp <- tryCatch(fread(file), error = function(e) NULL)
    if (is.null(dt_tmp)) return(NULL)

    dt_tmp <- dt_tmp[
      Group %chin% keep_groups &
      max_MAF != 0.01 &
      Region %chin% regions
    ]

    if (nrow(dt_tmp) == 0L) return(NULL)

    dt_tmp[, `:=`(
      phenotype = pheno,
      ancestry = ancestry,
      dataset = basename(info$dataset)
    )]

    dt_tmp
  }, future.seed = TRUE)
})

dt_biobank <- rbindlist(dt_list, fill = TRUE, use.names = TRUE)
fwrite(dt_biobank, sep = "\t", quote=FALSE,
  file = "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/biobank_results_at_sig_genes.tsv.gz")

# Finally, for the table in the paper of significant results - let's merge in the manually curated data as well.
# First, download the latest version of the gene-phenotype pairs
# (this code must to run locally as we are downloading from google sheets)
system("scp qen698@cluster2.bmrc.ox.ac.uk:/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent_summary.tsv.gz manuscript_tables/Tables/")
dt <- fread("manuscript_tables/Tables/gene_phenotype_pairs_for_agent_summary.tsv.gz")

# To do this we will grab the existing version and add merge to update the
# collection, and spit out a tsv to be re-uploaded.
library(googlesheets4)
# Read in the continuous manually curated traits
brava_supp_tables_url <- "https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=1082845169#gid=1082845169"
dt_cts <- read_sheet(brava_supp_tables_url, sheet="Table S12")
dt_binary <- read_sheet(brava_supp_tables_url, sheet="Table S13")
dt_cts <- data.table(dt_cts %>% 
  dplyr::select(`Phenotype ID`, `Gene ID`, `Significant in Genebass`,
    `Gene-bass significant pLoF`, `Gene-bass significant missense`))
dt_binary <- data.table(dt_binary %>% 
  dplyr::select(`Phenotype ID`, `Gene ID`,
    `Significant in Genebass`, `Gene-bass significant pLoF`,
    `Gene-bass significant missense`, `Gene OMIM (MIM number)`,
    `OMIM phenotype match`,
    `OMIM (MIM number) of phenotypes associated to genetic variation in the gene`,
    `Disease name of phenotypes associated to variation in the gene in OMIM`))

# Read in the binary manually curated traits
setkeyv(dt_cts, c("Phenotype ID", "Gene ID"))
setkeyv(dt_binary, c("Phenotype ID", "Gene ID"))
setkeyv(dt, c("Phenotype ID", "Gene ID"))

dt_cts_updated <- merge(dt %>% dplyr::filter(`Binary or continuous` == "continuous"), dt_cts, all.x=TRUE) %>% 
  dplyr::select(-`Binary or continuous`)
setcolorder(dt_cts_updated,
  c("Description", "Gene symbol", "Phenotype ID", "Gene ID", "Phenotype broad category",
    "N significant Cauchy associations across meta-analyses",
    "Significant Cauchy", "Min P-value of significant Cauchy associations",
    "Cauchy combination with minimum P-value",
    "N significant associations across (max MAF, mask, meta-analyses)",
    "N significant meta-analyses", "Significant meta-analyses", "Significant burden meta-analyses",
    "Significant SKAT/SKAT-O meta-analyses", "Significant masks", "Significant max MAFs",
    "Significant classes", "Meta-analysis method", "Min P-value of significant associations",
    "Meta-analysis with minimum P-value", "Mask for most significant association",
    "Any burden significant", "Any SKAT/SKAT-O significant", "Significant in Genebass",
    "Gene-bass significant pLoF", "Gene-bass significant missense", "Literature verdict",
    "Open targets verdict", "Max verdict", "Justification" ,"Citations (PubMed IDs)"))

# Ensure that all of the required columns are in there - and then fill in the remainder, similar for the binary traits.
# Merge in any new results gathered from the novelty agent (if it was required to be run for new associations)
fwrite(dt_cts_updated, "cts_supp_table.tsv", sep='\t')

dt_binary_updated <- merge(dt %>% dplyr::filter(`Binary or continuous` == "binary"), dt_binary, all.x=TRUE) %>%
  dplyr::select(-`Binary or continuous`)
setcolorder(dt_binary_updated,
  c("Description", "Gene symbol", "Phenotype ID", "Gene ID", "Phenotype broad category",
    "N significant Cauchy associations across meta-analyses",
    "Significant Cauchy", "Min P-value of significant Cauchy associations",
    "Cauchy combination with minimum P-value",
    "N significant associations across (max MAF, mask, meta-analyses)",
    "N significant meta-analyses", "Significant meta-analyses", "Significant burden meta-analyses",
    "Significant SKAT/SKAT-O meta-analyses", "Significant masks", "Significant max MAFs",
    "Significant classes", "Meta-analysis method", "Min P-value of significant associations",
    "Meta-analysis with minimum P-value", "Mask for most significant association",
    "Any burden significant", "Any SKAT/SKAT-O significant", "Significant in Genebass",
    "Gene-bass significant pLoF", "Gene-bass significant missense",
    "Gene OMIM (MIM number)", "OMIM phenotype match",
    "OMIM (MIM number) of phenotypes associated to genetic variation in the gene",
    "Disease name of phenotypes associated to variation in the gene in OMIM", 
    "Literature verdict", "Open targets verdict", "Max verdict", "Justification",
    "Citations (PubMed IDs)"))
# Merge in any new novelty agent results that we need
fwrite(dt_binary_updated, "binary_supp_table.tsv", sep='\t')

# This can now be pasted back in to the supplementary table google sheets and the final 
# manual pieces of information filled in.

# For Tables S14 and S15, we can just take the above before the summary part (keep it long)
dt <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs_for_agent.tsv.gz") %>%
  rename(Description = phenotype_full,
      `Gene symbol` = external_gene_name,
      `Phenotype ID` = phenotype,
      `Gene ID` = Region,
      `Phenotype broad category` = phenotype_broad_category,
      `meta analyzed` = meta,
      `Mask` = Group,
      `max MAF` = max_MAF,
      class = class,
      `Meta analysis method` = type,
      Pvalue = Pvalue, 
      `BETA Burden` = BETA_Burden, 
      `SE Burden` = SE_Burden)
setcolorder(dt, c("Description", "Gene symbol", "Phenotype ID",
  "Gene ID", "Phenotype broad category", "meta analyzed",
  "Mask", "max MAF", "class", "Meta analysis method",
  "Pvalue", "BETA Burden", "SE Burden", "phenotype_class"))

fwrite(dt %>% filter(phenotype_class == "continuous") %>% select(-phenotype_class),
  file = "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/cts_long.tsv", sep='\t', quote=FALSE)
   
fwrite(dt %>% filter(phenotype_class == "binary") %>% select(-phenotype_class),
  file = "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/binary_long.tsv", sep='\t', quote=FALSE)

# for (file in files) {
#   cat(file, "...\n")
#   binary <- grepl("binary", file)
#   dt_filtered <- fread(file)
#   dt_all <- fread("manuscript_figures/biobank_results_at_sig_genes.tsv.gz") %>% 
#     filter(Pvalue_SKAT < 2.5e-7 | Pvalue < 2.5e-7 | Pvalue_Burden < 6.7e-7)
#   dt_meta <- fread("manuscript_figures/gene_phenotype_pairs_for_agent_superset.tsv.gz") %>% 
#     filter(((class %in% c("SKAT-O", "SKAT") & Pvalue < 2.5e-7)) |
#       ((class == "Burden") & Pvalue < 6.7e-7)) %>% 
#     select(Region, Group, max_MAF, class, phenotype, meta)
#   dt_all_long <- melt(
#     dt_all,
#     measure.vars = c("Pvalue", "Pvalue_Burden", "Pvalue_SKAT"),
#     variable.name = "class",
#     value.name = "Pvalue"
#   )

#   # Map column names to desired labels
#   dt_all_long[, class := fifelse(
#     class == "Pvalue", "SKAT-O",
#     fifelse(class == "Pvalue_Burden", "Burden", "SKAT")
#   )]

#   dt_all_long_filter <- dt_all_long %>% select(Region, Group, max_MAF, class, phenotype)

#   meta_filter <- setdiff(dt_meta %>% select(-meta), dt_all_long_filter)
#   setkeyv(meta_filter, c("Region", "Group", "max_MAF", "class", "phenotype"))
#   meta_filter <- merge(dt_meta, meta_filter) %>% 
#     mutate(test = paste0(
#       "meta=", meta, ":",
#       class, ":",
#       ifelse(class == "Burden", "Inverse variance weighted", "Stouffer"), ":",
#       Group, ":",
#       "maf=", max_MAF)) %>% group_by(Region, phenotype) %>% summarise(tests = paste(test, collapse=", ")) %>%
#     rename(`Gene ID`=Region, `Phenotype ID`=phenotype)
#   meta_filter <- data.table(meta_filter)

#   setkeyv(meta_filter, c("Gene ID", "Phenotype ID"))
#   setkeyv(dt_filtered, c("Gene ID", "Phenotype ID"))
#   # Then merge in the results from dt_binary and dt_cts
#   dt_final <- merge(meta_filter, dt_filtered)
#   dt_final_filtered <- dt_final %>% filter(`Max verdict` %in% c("Not found", "Hypothesized"))
#   print(dt_final_filtered)
#   # fwrite(dt_final_filtered, file=gsub(".tsv", "_final_filtered.tsv", file), sep='\t')
#   # fwrite(dt_final_filtered %>% filter(grepl("=ALL", tests)), file=gsub(".tsv", "_final_filtered_ALL.tsv", file), sep='\t')
# }





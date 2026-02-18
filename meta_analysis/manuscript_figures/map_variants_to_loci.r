#!/usr/bin/env Rscript
# map_variants_to_loci.R
#
# Usage:
# Rscript map_variants_to_loci.R --loci loci.tsv[.gz] --variants variants.tsv[.gz] --out mapped.tsv --byPhen TRUE
#
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

option_list <- list(
  make_option(c("--loci"), type="character", help="Loci TSV (must contain locus_id, locus_start, locus_end, chr)"),
  make_option(c("--variants"), type="character", help="Variants TSV (must contain #CHROM, POS; optional REF, ALT, phenotypeID)"),
  make_option(c("--out"), type="character", default="mapped_variants_to_loci.tsv",
              help="Output TSV filename"),
  make_option(c("--byPhen"), type="logical", default=FALSE,
              help="TRUE: only match variants to loci with same phenotypeID (requires 'phenotype' column in both files). FALSE: match to any locus"),
  make_option(c("--chr.col"), type="character", default="chr", help="Chromosome column name in loci file")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$loci) || is.null(opt$variants)) stop("Please supply --loci and --variants")

# helper: read with fread (handles gz)
read_dt <- function(path) {
  dt <- fread(path)
  if (nrow(dt) == 0) warning("Empty table: ", path)
  return(dt)
}

loci_dt <- read_dt(opt$loci)
vars_dt <- read_dt(opt$variants)

# minimal column checks
needed_loci <- c("locus_id","locus_start","locus_end", opt$`chr.col`)
if (!all(needed_loci %in% names(loci_dt))) {
  stop("Loci file must contain columns: ", paste(needed_loci, collapse=", "))
}
if (!("#CHROM" %in% names(vars_dt) && "POS" %in% names(vars_dt))) {
  stop("Variant file must contain columns: #CHROM  and POS")
}

# normalize chrom strings (remove leading 'chr' for matching consistency)
norm_chr <- function(x) sub("^chr", "", as.character(x), ignore.case = TRUE)
loci_dt[, chr_norm := norm_chr(get(opt$`chr.col`))]
vars_dt[, chr_norm := norm_chr(get("#CHROM"))]

# prepare loci as intervals
loci_dt[, locus_start := as.integer(locus_start)]
loci_dt[, locus_end := as.integer(locus_end)]
# ensure start <= end
loci_dt <- loci_dt[locus_end >= locus_start]

# prepare variants as intervals (point intervals pos..pos)
vars_dt[, pos_start := as.integer(get("POS"))]
vars_dt[, pos_end := pos_start]

setnames(loci_dt,
         old=c("chr_norm","locus_start","locus_end"),
         new=c("seqname","locus_start","locus_end"))

setnames(vars_dt,
         old=c("chr_norm","pos_start","pos_end"),
         new=c("seqname","var_start","var_end"))

# create copies to avoid clobbering original names
loci_intervals <- copy(loci_dt)
vars_intervals <- copy(vars_dt)

# Set keys for foverlaps
setkeyv(loci_dt, c("seqname", "locus_start", "locus_end"))
setkeyv(vars_dt, c("seqname", "var_start", "var_end"))

# if matching by phenotype, ensure column exists
if (opt$byPhen) {
  if (!("phenotype" %in% names(loci_intervals))) stop("Loci file must have 'phenotype' column for --byPhen=TRUE")
  if (!("phenotypeID" %in% names(vars_intervals))) stop("Variants file must have 'phenotypeID' column for --byPhen=TRUE")
  # we'll do phenotype-wise joins to be efficient
  message("Matching by phenotype; will join per-phenotype.")
  out_list <- list()
  phen_list <- intersect(unique(loci_intervals$phenotype), unique(vars_intervals$phenotypeID))
  if (length(phen_list) == 0) {
    warning("No overlapping phenotypes between loci and variants.")
  }
  for (ph in phen_list) {
    sub_l <- loci_intervals[phenotype == ph]
    sub_v <- vars_intervals[phenotypeID == ph]
    if (nrow(sub_l) == 0 || nrow(sub_v) == 0) next
    res <- foverlaps(sub_v, sub_l,
      by.x = c("seqname","var_start","var_end"),
      by.y = c("seqname","locus_start","locus_end"),
      nomatch = 0L)
    out_list[[ph]] <- res
  }
  if (length(out_list) == 0) {
    mapped <- data.table() 
  } else {
    mapped <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
  }
} else {
  # match to any locus
  message("Matching variants to any locus (no phenotype constraint).")
  mapped <- foverlaps(vars_intervals, loci_intervals,
    by.x = c("seqname","var_start","var_end"),
    by.y = c("seqname","locus_start","locus_end"),
    nomatch = 0L)
}

# Post-process: restore original column names and tidy output
if (nrow(mapped) == 0) {
  message("No variants mapped to loci based on the provided files/criteria.")
  fwrite(data.table(), file = opt$out, sep = "\t")
} else {
  fwrite(out_dt, file = opt$out, sep = "\t", na = "NA", quote = FALSE)
  message("Wrote ", nrow(out_dt), " mapped rows to ", opt$out)
}

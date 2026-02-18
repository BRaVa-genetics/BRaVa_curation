#!/usr/bin/env Rscript
#
# count_phen_loci.R
#
# Usage:
#   Rscript count_phen_loci.R --in variants.tsv --out mapped_loci.tsv --window 500000
#
# Input must contain columns: phenotype, chr, pos
# (optional: ref, alt)
#
# Output: prints per-phenotype locus counts and writes mapping file if --out provided.
#

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  # GenomicRanges is Bioconductor; attempt to load
  have_gr <- requireNamespace("GenomicRanges", quietly = TRUE)
  if (!have_gr) {
    message("Package 'GenomicRanges' not found. The script will try to install it now (requires BiocManager).")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("GenomicRanges", ask = FALSE)
    have_gr <- requireNamespace("GenomicRanges", quietly = TRUE)
    if (!have_gr) stop("Failed to install/load GenomicRanges. Please install it manually and rerun.")
  }
  library(GenomicRanges)
})

option_list <- list(
  make_option(c("-i","--infile"), type="character", help="Input TSV/CSV file with header (phenotypeID, #CHROM, POS, optional REF,ALT)"),
  make_option(c("-o","--out"), type="character", default=NULL, help="Optional output TSV mapping variants->loci"),
  make_option(c("-w","--window"), type="integer", default=500000, help="Window size in bp (default 500000)"),
  make_option(c("-s","--sep"), type="character", default=NULL, help="Separator for input file (auto-detected if not provided)")
)

opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$infile)) stop("Please provide --in <variants.tsv>")

# read input
dt <- fread(opt$infile, sep = opt$sep)
required_cols <- c("phenotypeID","#CHROM","POS")
if (!all(required_cols %in% names(dt))) {
  stop("Input must contain columns: phenotype, chr, pos. Found: ", paste(names(dt), collapse=", "))
}

# normalize chromosome string (remove leading 'chr' if present)
dt[, chr2 := sub("^chr", "", `#CHROM`, ignore.case = TRUE)]
# ensure pos numeric
dt[, pos := as.integer(POS)]

# We'll process per phenotype -> create GRanges with ranges +/- window
window <- as.integer(opt$window)

# Prepare result container
mapping_list <- list()
summary_list <- list()
locus_global_counter <- 0L

phenotypes <- unique(dt$phenotypeID)
setorder(dt, phenotypeID, chr2, POS)

for (phen in phenotypes) {
  dtp <- dt[phenotypeID == phen]
  # process per chromosome
  loci_for_phen <- list()
  for (chrom in unique(dtp$chr2)) {
    sub <- dtp[chr2 == chrom]
    starts <- pmax(1L, sub$POS - window)
    ends   <- sub$POS + window
    gr <- GenomicRanges::GRanges(seqnames = chrom,
                                 ranges = IRanges::IRanges(start = starts, end = ends))
    # reduce to merge overlapping intervals (loci)
    gr_reduced <- GenomicRanges::reduce(gr, ignore.strand = TRUE)
    if (length(gr_reduced) == 0) next
    # assign locus ids for this reduced set
    # find overlaps of original variants to reduced intervals
    hits <- GenomicRanges::findOverlaps(gr, gr_reduced, ignore.strand = TRUE)
    # hits maps original interval index -> reduced index
    # build mapping table
    df_hits <- data.table(orig_idx = queryHits(hits),
                          reduced_idx = subjectHits(hits))
    # for each reduced interval create a locus id and record members
    for (ridx in seq_len(length(gr_reduced))) {
      member_rows <- df_hits[reduced_idx == ridx, orig_idx]
      if (length(member_rows) == 0) next
      locus_global_counter <- locus_global_counter + 1L
      locus_id <- paste0("L", locus_global_counter)
      locus_start <- start(gr_reduced[ridx])
      locus_end   <- end(gr_reduced[ridx])
      # original row indices in 'sub' are 1..nrow(sub) corresponding to orig_idx
      members_dt <- data.table(
        phenotype = phen,
        chr = chrom,
        pos = sub$pos[member_rows],
        ref = if ("ref" %in% names(sub)) sub$ref[member_rows] else NA_character_,
        alt = if ("alt" %in% names(sub)) sub$alt[member_rows] else NA_character_,
        locus_id = locus_id,
        locus_start = locus_start,
        locus_end = locus_end
      )
      mapping_list[[length(mapping_list)+1]] <- members_dt
    }
  }

  # count loci for this phenotype: derived from how many distinct locus_ids we assigned for this phenotype
  if (length(mapping_list) > 0) {
    # count how many locus_ids were created for this phenotype by checking mapping_list items
    # but easier: count locus ids assigned so far minus previous counter snapshot
    # We'll compute per-phenotype counts by inspecting mapping_list entries mode:
    # collect locus_ids created for this phenotype:
    all_map_for_phen <- rbindlist(mapping_list, fill=TRUE)
    n_loci_phen <- length(unique(all_map_for_phen[phenotype == phen, locus_id]))
  } else {
    n_loci_phen <- 0L
  }
  summary_list[[as.character(phen)]] <- n_loci_phen
}

# Consolidate mapping and summary
if (length(mapping_list) > 0) {
  mapping_dt <- rbindlist(mapping_list, fill=TRUE)
  # reorder columns
  setcolorder(mapping_dt, c("phenotype","chr","pos","ref","alt","locus_id","locus_start","locus_end"))
} else {
  mapping_dt <- data.table(phenotype=character(), chr=character(), pos=integer(),
                           ref=character(), alt=character(), locus_id=character(),
                           locus_start=integer(), locus_end=integer())
}

# compute per-phenotype counts robustly from mapping_dt
phen_counts <- mapping_dt[, .(n_loci = uniqueN(locus_id)), by = phenotype]
setorder(phen_counts, -n_loci)

# Print summary
cat("Per-phenotype locus counts (locus = variant +/-", window, "bp):\n")
print(phen_counts)
total_pairs <- sum(phen_counts$n_loci)
cat("\nTotal phenotype-locus pairs:", total_pairs, "\n")

# write mapping file if requested
if (!is.null(opt$out)) {
  fwrite(mapping_dt, file = opt$out, sep = "\t", na = "NA", quote = FALSE)
  cat("Wrote mapping of variants -> loci to", opt$out, "\n")
}

#!/usr/bin/env Rscript
library(data.table)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(ggtext)
library(ggnewscale)

source("meta_analysis_utils.r")

local <- TRUE
if (local) {
  data_dir <- "manuscript_figures/"
  out_dir <- "manuscript_figures/Figures/"
} else {
  data_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/"
  out_dir <- "/well/lindgren/dpalmer/BRaVa_curation/meta_analysis/manuscript_figures/Figures/"
}

thres <- 0.05 / (20000 * 3 * 2)

dt <- fread(paste0(data_dir, "biobank_results_at_sig_genes.tsv.gz"))
setDT(dt)

dt_agent <- fread(paste0(data_dir, "novelty_agent_results.tsv"))
# Ensure that the chronic kidney disease is correctly named
dt_agent <- dt_agent %>% mutate(phenotype = ifelse(phenotype == "Chronic renal failure", "Chronic kidney disease", phenotype))

setDT(dt_agent)
setnames(dt_agent, old = c("gene", "phenotype"), new = c("external_gene_name", "phenotype_full"))
dt_agent <- dt_agent[!(external_gene_name %in% c("DNMT3A", "TET2", "ASXL1"))]
dt_agent[, pair := paste0(phenotype_full, " - ", external_gene_name)]
bottom_df_all <- unique(dt_agent[, .(pair, max_verdict)])

dt_meta_superset_all <- fread(paste0(data_dir, "gene_phenotype_pairs_for_agent_superset.tsv.gz"))
dt_meta_superset_all <- dt_meta_superset_all %>% mutate(phenotype_full = unlist(renaming_phenotype_list[phenotype]))
setDT(dt_meta_superset_all)

# dt_meta_subset_all <- fread(paste0(data_dir, "gene_phenotype_pairs_for_agent.tsv.gz"))
dt_meta_subset_all <- fread(paste0(data_dir, "gene_phenotype_pairs_for_agent_cauchy.tsv.gz"))
dt_meta_subset_all <- dt_meta_subset_all %>% mutate(phenotype_full = unlist(renaming_phenotype_list[phenotype]))
dt_meta_subset_all <- dt_meta_subset_all[(class == "Burden" & Pvalue < thres)]

setDT(dt_meta_subset_all)
dt_meta_subset_all <- dt_meta_subset_all[!grepl("just_uk-biobank_and_all-of-us", file)]
dt_meta_subset_all[external_gene_name == "", external_gene_name := Region]

keep_pairs <- dt[, .(any_sig = any(Pvalue_Burden < thres, na.rm = TRUE)),
                 by = .(Region, phenotype)][any_sig == FALSE, .(Region, phenotype)]
setkey(keep_pairs, Region, phenotype)

color_amr <- "#ED1E24"
color_eur <- "#6AA5CD"
color_afr <- "#941494"
color_sas <- "#FF9912"
color_eas <- "#108C44"
color_oth <- "#ABB9B9"
color_mde <- "#33CC33"
color_asj <- "coral"
color_nfe <- color_eur
color_fin <- "#002F6C"

pop_colors <- c(
  AFR = color_afr,
  AMR = color_amr,
  EAS = color_eas,
  FIN = color_fin,
  EUR = color_nfe,
  NEF = color_nfe,
  OTH = color_oth,
  SAS = color_sas,
  MDE = color_mde,
  ASJ = color_asj,
  uniform = "pink",
  consanguineous = "pink",
  SAS_non_consang = "orange",
  ALL = "black",
  `non-EUR` = color_oth
)

verdict_cols <- c(
  "Not found" = "#D73027",
  "Hypothesized" = "#FC8D59",
  "Existing" = "#FEE090",
  "Established" = "#FFFFBF",
  "Unknown" = "grey85"
)

fill_vals <- c(
  "Significant (< 4.2 × 10⁻⁷)" = "#084081",
  "< 1 × 10⁻⁴" = "#2b8cbe",
  "< 0.05" = "#c7eae5",
  "> 0.05" = "grey90"
)

ancestry_order <- c("AFR", "AMR", "EAS", "EUR", "SAS", "non-EUR", "ALL")
fixed_point_size <- 4.0

style_one_label <- function(plain_label) {
  parts <- strsplit(plain_label, ",\\s*", perl = TRUE)[[1]]
  dataset <- parts[1]
  ancestry <- if (length(parts) >= 2) parts[length(parts)] else ""

  anc_col <- "black"
  if (!is.na(ancestry) && ancestry != "" && ancestry %in% names(pop_colors)) {
    anc_col <- pop_colors[[ancestry]]
  }

  if (ancestry != "") {
    ancestry_html <- paste0(
      "<span style=\"color:", anc_col, ";\">",
      ancestry,
      "</span>"
    )
    label <- paste0(dataset, ", ", ancestry_html)
  } else {
    label <- dataset
  }

  if (grepl("(?i)Meta-analysis", dataset, perl = TRUE)) {
    label <- paste0("<b>", label, "</b>")
  }

  if (plain_label == "Novelty / Agent") {
    label <- paste0("<i>", label, "</i>")
  }

  label
}

library(future.apply)
library(progressr)

plan(multisession, workers = 4)   # adjust as needed
handlers("progress")

project_dir <- getwd()

# Base jobs: one set of temp files per type × maf × G
base_job_dt <- CJ(
  type = c("binary", "continuous"),
  maf = c(0.001, 0.0001),
  G = c(
    "pLoF",
    "damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering"
  )
)

tmp_dir <- file.path(tempdir(), "gene_pheno_plot_jobs")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

base_job_dt[, `:=`(
  dt_file = file.path(
    tmp_dir,
    paste0(
      "dt_",
      type, "_",
      maf, "_",
      gsub("[^A-Za-z0-9]+", "_", G),
      ".tsv.gz"
    )
  ),
  dt_meta_file = file.path(
    tmp_dir,
    paste0(
      "dt_meta_",
      type, "_",
      maf, "_",
      gsub("[^A-Za-z0-9]+", "_", G),
      ".tsv.gz"
    )
  ),
  dt_merge_file = file.path(
    tmp_dir,
    paste0(
      "dt_merge_",
      type, "_",
      maf, "_",
      gsub("[^A-Za-z0-9]+", "_", G),
      ".tsv.gz"
    )
  )
)]

# Write each base subset once in the main session
for (i in seq_len(nrow(base_job_dt))) {
  type_i <- base_job_dt$type[i]
  maf_i  <- base_job_dt$maf[i]
  G_i    <- base_job_dt$G[i]

  dt_sub <- dt[max_MAF == maf_i & Group == G_i]
  fwrite(dt_sub, base_job_dt$dt_file[i])

  dt_meta_sub <- dt_meta_superset_all[max_MAF == maf_i & Group == G_i & class == "Burden"]
  drop_cols <- intersect(
    names(dt_meta_sub),
    c("external_gene_name", "phenotype_full", "phenotype_broad_category", "phenotype_class")
  )
  if (length(drop_cols)) dt_meta_sub[, (drop_cols) := NULL]
  fwrite(dt_meta_sub, base_job_dt$dt_meta_file[i])

  dt_merge_in_sub <- unique(
    dt_meta_subset_all[
      max_MAF == maf_i & Group == G_i & class == "Burden",
      .(Region, external_gene_name, phenotype_full, phenotype, phenotype_class, phenotype_broad_category)
    ]
  )
  fwrite(dt_merge_in_sub, base_job_dt$dt_merge_file[i])
}

# Add include_height as a separate job dimension
job_dt <- CJ(
  type = c("binary", "continuous"),
  maf = c(0.001, 0.0001),
  G = c(
    "pLoF",
    "damaging_missense_or_protein_altering",
    "pLoF;damaging_missense_or_protein_altering"
  ),
  include_height = c(TRUE, FALSE)
)

# Remove impossible combination: binary + no height
job_dt <- job_dt[!(type == "binary" & include_height == FALSE)]

# Merge in the file paths from the base jobs
setkey(base_job_dt, type, maf, G)
setkey(job_dt, type, maf, G)
job_dt <- base_job_dt[job_dt, on = .(type, maf, G)]

with_progress({
  p <- progressor(along = seq_len(nrow(job_dt)))

  future_lapply(seq_len(nrow(job_dt)), function(i) {

    setwd(project_dir)

    type_i <- job_dt$type[i]
    maf_i  <- job_dt$maf[i]
    G_i    <- job_dt$G[i]
    include_height_i <- job_dt$include_height[i]

    p(sprintf(
      "%s | maf=%s | %s | %s",
      type_i,
      maf_i,
      G_i,
      ifelse(include_height_i, "with_height", "no_height")
    ))

    dt_tmp <- fread(job_dt$dt_file[i])
    dt_meta <- fread(job_dt$dt_meta_file[i])
    dt_merge_in <- fread(job_dt$dt_merge_file[i])

    if (!include_height_i) {
      dt_tmp <- dt_tmp[phenotype != "Height"]
      dt_meta <- dt_meta[phenotype != "Height"]
      dt_merge_in <- dt_merge_in[phenotype != "Height"]
    }

    setkey(dt_merge_in, phenotype, Region)
    setkey(dt_meta, phenotype, Region)
    dt_meta <- merge(dt_meta, dt_merge_in, all.x = TRUE, sort = FALSE)

    setkey(dt_tmp, phenotype, Region)
    dt_tmp <- merge(dt_tmp, dt_merge_in, all.x = TRUE, sort = FALSE)
    dt_tmp[external_gene_name == "", external_gene_name := Region]

    dt_tmp <- dt_tmp[keep_pairs, on = .(Region, phenotype), nomatch = 0L]

    dt_tmp[, Pvalue := Pvalue_Burden]
    drop_cols <- intersect(names(dt_tmp), c("Pvalue_Burden", "Pvalue_SKAT"))
    if (length(drop_cols)) dt_tmp[, (drop_cols) := NULL]
    dt_tmp[, dataset := unlist(renaming_plot_biobank_list[dataset])]

    dt_meta[, `:=`(
      dataset = "Meta-analysis",
      ancestry = meta
    )]
    dt_meta <- dt_meta[keep_pairs, on = .(Region, phenotype), nomatch = 0L]

    common_cols <- intersect(names(dt_tmp), names(dt_meta))
    dt_tmp <- rbindlist(
      list(dt_tmp[, ..common_cols], dt_meta[, ..common_cols]),
      use.names = TRUE,
      fill = TRUE
    )

    dt_tmp[, pair := paste0(phenotype_full, " - ", external_gene_name)]

    plot_dt <- dt_tmp[phenotype_class == type_i]
    plot_dt[, ancestry := fifelse(ancestry == "non_EUR", "non-EUR", ancestry)]
    plot_dt[, row_label := paste0(dataset, ", ", ancestry)]
    plot_dt[, row_label := factor(row_label)]
    plot_dt[, pair := factor(pair)]

    plot_dt[, p_bin := factor(
      cut(
        Pvalue,
        breaks = c(-Inf, thres, 1e-4, 0.05, Inf),
        labels = c(
          "Significant (< 4.2 × 10⁻⁷)",
          "< 1 × 10⁻⁴",
          "< 0.05",
          "> 0.05"
        ),
        right = FALSE
      ),
      levels = names(fill_vals)
    )]

    plot_dt[, sign_char := fifelse(
      is.na(BETA_Burden), "",
      fifelse(BETA_Burden > 0, "+",
              fifelse(BETA_Burden < 0, "\u2212", ""))
    )]

    dark_bins <- c("Significant (< 4.2 × 10⁻⁷)", "< 1 × 10⁻⁴")
    plot_dt[, txt_col := fifelse(p_bin %in% dark_bins, "white", "black")]

    row_levels <- unique(
      plot_dt[, .(dataset, ancestry, row_label)][
        , ancestry := factor(ancestry, levels = ancestry_order)
      ][
        order(ancestry, dataset)
      ]$row_label
    )
    row_levels <- as.character(c(
      row_levels[!grepl("Meta-analysis", row_levels)],
      row_levels[grepl("Meta-analysis", row_levels)]
    ))

    bottom_df <- unique(bottom_df_all[pair %in% levels(plot_dt$pair)])
    # Remove pairs with no novelty annotation
    bottom_df <- bottom_df[!is.na(max_verdict)]

    # Keep only the plot rows that still have novelty-agent info
    valid_pairs <- unique(bottom_df$pair)
    plot_dt <- plot_dt[pair %in% valid_pairs]
    plot_dt[, pair := factor(pair)]

    # Rebuild bottom_df so its factor levels match the filtered plot_dt
    bottom_df <- bottom_df[pair %in% valid_pairs]
    bottom_df[, pair := factor(pair, levels = levels(plot_dt$pair))]

    bottom_row_label <- "Novelty / Agent"
    new_levels <- c(row_levels, bottom_row_label)
    plot_dt[, row_label := factor(row_label, levels = rev(new_levels))]
    bottom_df[, row := factor(bottom_row_label, levels = levels(plot_dt$row_label))]

    label_spacers <- c("LABEL1", "LABEL2", "LABEL3", "LABEL4")
    x_levels <- c(label_spacers, levels(plot_dt$pair))

    plot_dt[, pair := factor(as.character(pair), levels = x_levels)]
    bottom_df[, pair := factor(as.character(pair), levels = x_levels)]

    label_df <- data.table(
      row_label = factor(levels(plot_dt$row_label), levels = levels(plot_dt$row_label)),
      pair = factor("LABEL4", levels = x_levels),
      label_html = vapply(levels(plot_dt$row_label), style_one_label, FUN.VALUE = character(1))
    )

    p <- ggplot() +
      geom_tile(
        data = bottom_df,
        aes(x = pair, y = row, fill = max_verdict),
        width = 0.95, height = 0.95
      ) +
      scale_fill_manual(name = "Max verdict", values = verdict_cols, na.value = "grey85") +
      ggnewscale::new_scale_fill() +
      geom_point(
        data = plot_dt,
        aes(x = pair, y = row_label, fill = p_bin),
        shape = 21, color = "black", stroke = 0.25, size = fixed_point_size
      ) +
      geom_text(
        data = plot_dt,
        aes(x = pair, y = row_label, label = sign_char, color = txt_col),
        size = 3, fontface = "bold", vjust = 0.5
      ) +
      geom_richtext(
        data = label_df,
        aes(x = pair, y = row_label, label = label_html),
        inherit.aes = FALSE,
        hjust = 1,
        vjust = 0.5,
        fill = NA,
        label.color = NA,
        label.padding = grid::unit(c(0, 0, 0, 0), "pt"),
        label.margin = grid::unit(c(0, 0, 0, 0), "pt"),
        label.r = grid::unit(0, "pt"),
        size = 2.8
      ) +
      scale_fill_manual(
        values = fill_vals,
        name = expression("Burden " * italic(P) * "-value"),
        na.value = "grey90",
        drop = FALSE,
        limits = names(fill_vals)
      ) +
      scale_color_identity() +
      scale_x_discrete(
        limits = x_levels,
        labels = function(x) {
          labs <- vector("list", length(x))

          for (j in seq_along(x)) {
            lbl <- x[j]

            if (lbl %in% label_spacers) {
              labs[[j]] <- quote("")
              next
            }

            parts <- strsplit(lbl, " - ", fixed = TRUE)[[1]]

            if (length(parts) < 2) {
              labs[[j]] <- bquote(.(lbl))
              next
            }

            phenotype <- paste(parts[-length(parts)], collapse = " - ")
            gene <- parts[length(parts)]

            labs[[j]] <- bquote(.(phenotype) * ": " * italic(.(gene)))
          }

          as.expression(labs)
        },
        expand = expansion(add = c(0.03, 0.03)),
        drop = FALSE
      ) +
      scale_y_discrete(expand = expansion(add = c(0.03, 0.03))) +
      labs(
        x = NULL, y = NULL,
        title = ifelse(
          G_i == "damaging_missense_or_protein_altering",
          "Damaging missense or protein altering",
          ifelse(
            G_i == "pLoF;damaging_missense_or_protein_altering",
            "pLoF; Damaging missense or protein altering",
            G_i
          )
        ),
        subtitle = paste0("max MAF = ", ifelse(maf_i == 0.0001, "1 × 10⁻⁴", maf_i))
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title.position = "plot",
        plot.title = element_text(hjust = 0, face = "bold"),
        plot.subtitle = element_text(hjust = 0),
        plot.margin = grid::unit(c(1.8, 0.6, 0.6, 0.6), "cm")
      ) +
      guides(fill = guide_legend(order = 1)) +
      coord_fixed(ratio = 1, clip = "off")

    ncols <- length(x_levels)
    nrows <- length(levels(plot_dt$row_label))
    cell_width_cm <- 0.5
    cell_height_cm <- 0.5

    left_margin_cm <- 6.5
    right_margin_cm <- 5.5
    top_margin_cm <- 3.5
    bottom_margin_cm <- 1.5

    width_cm <- ncols * cell_width_cm + left_margin_cm + right_margin_cm
    height_cm <- nrows * cell_height_cm + top_margin_cm + bottom_margin_cm
    width_in <- width_cm / 2.54
    height_in <- height_cm / 2.54

    out_file <- paste0(
      out_dir,
      maf_i,
      G_i,
      type_i,
      ifelse(include_height_i, "", "_no_height")
    )

    cairo_pdf(filename = paste0(out_file, ".pdf"), width = width_in, height = height_in)
    print(p)
    dev.off()
  }, future.seed = TRUE)
})



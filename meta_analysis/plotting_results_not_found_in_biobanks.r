#!/usr/bin/env Rscript
library(data.table)
library(biomaRt)
# required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(ggtext)
library(ggnewscale)

source("meta_analysis_utils.r")
# system("scp 'qen698@cluster2.bmrc.ox.ac.uk:/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/gene_phenotype_pairs*.gz' manuscript_figures/")
# system("scp 'qen698@cluster2.bmrc.ox.ac.uk:/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/biobank_results_at_sig_genes.tsv.gz' manuscript_figures/")
local <- TRUE
if (local) {
  data_dir <- "manuscript_figures/"
  out_dir <- "manuscript_figures/Figures/"
} else {
  data_dir <- "/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/"
  out_dir <- "/well/lindgren/dpalmer/BRaVa_curation/meta_analysis/manuscript_figures/Figures/"
}

dt <- fread(paste0(data_dir, "biobank_results_at_sig_genes.tsv.gz"))
dt_agent <- fread(paste0(data_dir, "novelty_agent_results.tsv")) %>% 
  rename(external_gene_name = gene, phenotype_full = phenotype) %>% 
  filter(!(external_gene_name %in% c("DNMT3A", "TET2", "ASXL1")))


for (type in c("binary", "continuous")) {
  for (maf in c(0.001, 0.0001)) {
    for (G in c("pLoF", "damaging_missense_or_protein_altering", "pLoF;damaging_missense_or_protein_altering")) {
      dt_tmp <- dt %>% filter(max_MAF == !!maf, Group == !!G)
      dt_meta <- fread(paste0(data_dir, "gene_phenotype_pairs_for_agent_superset.tsv.gz")) %>% 
        filter(max_MAF == !!maf, Group == !!G, class == "Burden")
      dt_meta_subset <- fread(paste0(data_dir, "gene_phenotype_pairs_for_agent.tsv.gz")) %>% 
        filter(!grepl("just_uk-biobank_and_all-of-us", file)) %>% 
        filter(class == "Burden", Group == !!G, max_MAF == !!maf) %>%
        mutate(external_gene_name = ifelse(external_gene_name == "", Region, external_gene_name))
      dt_merge_in <- unique(
        dt_meta_subset %>% 
        dplyr::select(
          Region,
          external_gene_name,
          phenotype_full,
          phenotype,
          phenotype_class,
          phenotype_broad_category)
        )
      setkeyv(dt_merge_in, c("phenotype", "Region"))
      setkeyv(dt_meta, c("phenotype", "Region"))
      dt_meta <- merge(dt_meta, dt_merge_in)
      setkeyv(dt_tmp, c("phenotype", "Region"))
      # Merge in the gene-symbols and phenotype names
      dt_merge_in <- unique(dt_meta %>% 
        dplyr::select(
          Region,
          external_gene_name,
          phenotype_full,
          phenotype,
          phenotype_class,
          phenotype_broad_category))
      setkeyv(dt_merge_in, c("phenotype", "Region"))
      dt_tmp <- merge(dt_tmp, dt_merge_in)
      dt_tmp <- dt_tmp %>% mutate(external_gene_name = ifelse(external_gene_name == "", Region, external_gene_name))

      # find (Region, phenotype) pairs that have NO significant rows
      thresh <- 6.7e-7
      keep_pairs <- dt[, .(any_sig = any(Pvalue_Burden < thresh, na.rm = TRUE)),
                       by = .(Region, phenotype)][any_sig == FALSE, .(Region, phenotype)]

      # subset original table to only those pairs
      dt_filtered <- merge(dt_tmp, keep_pairs, by = c("Region", "phenotype"))
      dt_tmp <- dt_filtered %>% dplyr::select(-c("Pvalue", "Pvalue_SKAT")) %>% rename(Pvalue = Pvalue_Burden)
      dt_tmp[, dataset := unlist(renaming_plot_biobank_list[dataset])]
      dt_meta <- dt_meta %>% mutate(dataset = "Meta-analysis") %>% rename(ancestry = meta)
      dt_meta <- merge(dt_meta, keep_pairs, by = c("Region", "phenotype"))
      # Now, merge in the dt_meta

      dt_tmp <- rbind(dt_tmp %>% dplyr::select(intersect(names(dt_tmp), names(dt_meta))),
        dt_meta %>% dplyr::select(intersect(names(dt_tmp), names(dt_meta))))

      # prepare pair and ordering as before
      dt_tmp[, pair := paste0(phenotype_full, " - ", external_gene_name)]
      plot_dt <- dt_tmp %>% filter(phenotype_class == !!type)
      plot_dt[, ancestry := ifelse(ancestry == "non_EUR", "non-EUR", ancestry)]
      plot_dt[, row_label := paste0(dataset, ", ", ancestry)]
      plot_dt[, row_label := factor(row_label)]
      plot_dt[, pair := factor(pair)]

      # p-value bins (same bins as before)
      plot_dt[, p_bin := cut(Pvalue,
                             breaks = c(-Inf, 6.7e-7, 1e-4, 0.05, Inf),
                             labels = c(
                              "Significant (< 6.7 × 10⁻⁷)",
                              "< 1 × 10⁻⁴", "< 0.05", "> 0.05"),
                             right = FALSE)]

      fill_vals <- c(
        "Significant (< 6.7 × 10⁻⁷)" = "#084081",
        "< 1 × 10⁻⁴" = "#2b8cbe",
        "< 0.05" = "#c7eae5",
        "> 0.05" = "grey90"
      )

      # create sign label for each point: "+" for positive beta, "-" for negative, "" otherwise
      plot_dt[, sign_char := ifelse(
        is.na(BETA_Burden), "",
        ifelse(BETA_Burden > 0, "+",
               ifelse(BETA_Burden < 0, "\u2212", ""))
      )]

      # choose text color so symbol is readable on the fill
      dark_bins <- c("Significant (< 6.7 × 10⁻⁷)", "< 1 × 10⁻⁴")
      plot_dt[, txt_col := ifelse(p_bin %in% dark_bins, "white", "black")]

      ancestry_order <- c("AFR", "AMR", "EAS", "EUR", "SAS", "non-EUR", "ALL")

      # build ordered row table
      row_levels <- unique(
        plot_dt[
          , .(dataset, ancestry, row_label)
        ][
          , ancestry := factor(ancestry, levels = ancestry_order)
        ][
          order(ancestry, dataset)
        ]$row_label
      )
      row_levels <- as.character(c(row_levels[-grep("Meta-analysis", row_levels)], row_levels[grep("Meta-analysis", row_levels)]))
      # reverse so first row appears at the top of the plot
      plot_dt[, row_label := factor(row_label, levels = rev(row_levels))]

      # ---- ADD BOTTOM "Novelty / Agent" ROW ------------------------------------------------
      # construct bottom_df from dt_agent: pair mapping must match the 'pair' used above
      # pair format: "phenotype_full - external_gene_name"
      dt_agent[, pair := paste0(phenotype_full, " - ", external_gene_name)]
      bottom_df <- unique(dt_agent[, .(pair, max_verdict)])
      # keep only pairs present in the current plot
      bottom_df <- bottom_df[pair %in% levels(plot_dt$pair)]
      # map NA -> "Unknown"
      bottom_df[is.na(max_verdict), max_verdict := "Unknown"]
      # define verdict colours (adjust palette as you like)
      verdict_cols <- c(
        "Not found"    = "#D73027",
        "Hypothesized" = "#FC8D59",
        "Existing"     = "#FEE090",
        "Established"  = "#FFFFBF",
        "Unknown"      = "grey85"
      )
      # make factor levels consistent
      bottom_df[, max_verdict := factor(max_verdict, levels = names(verdict_cols))]

      # pick bottom row label and append it to the row factor levels
      bottom_row_label <- "Novelty / Agent"
      # append the bottom label to existing levels (so it will be rendered at the bottom)
      new_levels <- c(levels(plot_dt$row_label), bottom_row_label)
      plot_dt[, row_label := factor(row_label, levels = new_levels)]
      # add a 'row' column to bottom_df that matches the factor we just added
      bottom_df[, row := bottom_row_label]
      bottom_df[, row := factor(row, levels = levels(plot_dt$row_label))]
      # ------------------------------------------------------------------------------------

      # fixed point size
      fixed_point_size <- 4.0
      plot_dt <- plot_dt %>% filter(!(external_gene_name %in% c("DNMT3A", "TET2", "ASXL1")))

      ## ---- styled & coloured y-axis labels (drop in after you've set plot_dt$row_label) ----

      # your colour map (copy/paste)
      color_amr <- '#ED1E24'
      color_eur <- '#6AA5CD'
      color_afr <- '#941494'
      color_sas <- '#FF9912'
      color_eas <- '#108C44'
      color_oth <- '#ABB9B9'
      color_mde <- '#33CC33'
      color_asj <- 'coral'
      color_nfe <- color_eur
      color_fin <- '#002F6C'

      pop_colors <- c('AFR' = color_afr,
                      'AMR' = color_amr,
                      'EAS' = color_eas,
                      'FIN' = color_fin,
                      'EUR' = color_nfe,
                      'NEF' = color_nfe,
                      'OTH' = color_oth,
                      'SAS' = color_sas,
                      'MDE' = color_mde,
                      'ASJ' = color_asj,
                      'uniform' = 'pink',
                      'consanguineous' = 'pink',
                      'SAS_non_consang' = 'orange',
                      'ALL' = 'black',
                      'non-EUR' = color_oth)

      # helper: build styled label for a plain row label like "Meta-analysis, EUR" or "uk-biobank, AFR"
      style_one_label <- function(plain_label) {
        # split at the last comma-space to get dataset and ancestry (works with "dataset, ancestry")
        parts <- strsplit(plain_label, ",\\s*", perl=TRUE)[[1]]
        dataset <- parts[1]
        ancestry <- if (length(parts) >= 2) parts[length(parts)] else ""
        # get color for ancestry; fallback to black if missing
        if (ancestry != "") {
          anc_col <- pop_colors[[ancestry]]
        } else {
          anc_col <- "black"
        }
        # coloured ancestry span
        if (ancestry != "") {
          ancestry_html <- paste0("<span style='color:", anc_col, "'>", ancestry, "</span>")
          base_label <- paste0(dataset, ", ", ancestry_html)
        } else {
          base_label <- dataset
        }
        # bold meta rows
        if (grepl("(?i)Meta-analysis", dataset, perl=TRUE)) {
          base_label <- paste0("<b>", base_label, "</b>")
        }
        # italicize final novelty row if present (plain_label compared to exact string)
        if (plain_label == "Novelty / Agent") {
          base_label <- paste0("<i>", base_label, "</i>")
        }
        base_label
      }

      # build vector of styled labels in the same order as the current row_label factor levels
      orig_levels <- levels(plot_dt$row_label)   # original ordered labels
      styled_levels <- vapply(orig_levels, style_one_label, FUN.VALUE = character(1), USE.NAMES = FALSE)

      # now create the column in plot_dt and factor it to preserve ordering
      # map plain label -> styled label
      label_map <- setNames(styled_levels, orig_levels)
      plot_dt[, row_label_styled := label_map[as.character(row_label)] ]
      plot_dt[, row_label_styled := factor(row_label_styled, levels = styled_levels)]

      # ----- build plot: tiles (bottom) + dots, with clip = "off" and small axis expansion -----
      p <- ggplot() +
        geom_tile(
          data = bottom_df,
          aes(x = factor(pair, levels = levels(plot_dt$pair)), y = row, fill = max_verdict),
          width = 0.95, height = 0.95
        ) +
        scale_fill_manual(name = "Max verdict", values = verdict_cols, na.value = "grey85") +
        ggnewscale::new_scale_fill() +
        geom_point(
          data = plot_dt,
          aes(x = pair, y = row_label_styled, fill = p_bin),
          shape = 21, color = "black", stroke = 0.25, size = fixed_point_size
        ) +
        geom_text(
          data = plot_dt,
          aes(x = pair, y = row_label_styled, label = sign_char, color = txt_col),
          size = 3, fontface = "bold", vjust = 0.5
        ) +
        scale_fill_manual(values = fill_vals, name = "Burden P-value", na.value = "grey90", drop = FALSE) +
        scale_color_identity() +
        labs(x = NULL, y = NULL,
          title = ifelse(
            G == "damaging_missense_or_protein_altering", 
            "Damaging missense or protein altering",
            ifelse(
              G == "pLoF;damaging_missense_or_protein_altering",
              "pLoF; Damaging missense or protein altering",
              G)),
          subtitle = paste0("max MAF = ", ifelse(maf == 0.0001, "1 × 10⁻⁴", maf))) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
          axis.text.y = ggtext::element_markdown(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          # reserve extra top margin to be safe
          plot.margin = grid::unit(c(1.8, 0.6, 0.6, 0.6), "cm")
        ) +
        guides(fill = guide_legend(order = 1)) +
        # give a tiny bit of breathing room so edge cells aren't flush
        scale_x_discrete(expand = expansion(add = c(0.03, 0.03))) +
        scale_y_discrete(expand = expansion(add = c(0.03, 0.03))) +
        # allow geoms/text to draw beyond the panel limits (prevents clipping)
        coord_fixed(ratio = 1, clip = "off")

      # ----- compute figure size: increase top margin in paper size (cm) if still needed -----
      ncols <- length(levels(plot_dt$pair))
      nrows <- length(levels(plot_dt$row_label))
      cell_width_cm  <- 0.5
      cell_height_cm <- 0.5

      left_margin_cm  <- 2.5
      right_margin_cm <- 5.5
      top_margin_cm   <- 3.5   # increase further if still clipped
      bottom_margin_cm<- 1.5

      width_cm  <- ncols * cell_width_cm  + left_margin_cm + right_margin_cm
      height_cm <- nrows * cell_height_cm + top_margin_cm  + bottom_margin_cm
      width_in  <- (width_cm / 2.54)
      height_in <- (height_cm / 2.54)

      out_file <- paste0(out_dir, maf, G, type)
      cairo_pdf(filename = paste0(out_file, ".pdf"), width = width_in, height = height_in)
      print(p)
      dev.off()
      ggsave(plot=p, filename=paste0(out_file, ".png"), width = width_in, height = height_in)
    }
  }
}


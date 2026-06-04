#!/bin/Rscript
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)
library(eulerr)
library(colorspace)
library(ggtext)

source("../meta_analysis_utils.r")  # must include pop_colors here

# cauchy <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_cauchy_assocs_from_all_meta_subsets_101025_LOO.tsv.gz") %>% 
cauchy <- fread("../significant_cauchy_assocs_from_all_meta_subsets_101025_LOO.tsv.gz") %>% 
  filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456")))
cauchy$source <- "Cauchy"

cauchy <- cauchy %>% mutate(biobank = ifelse(biobank == "ALL", "ALL", unlist(renaming_plot_biobank_list[biobank])))

biobank_colors <- c(
  "All of Us"          = "#8CD1C5",
  "Biobank Japan"      = "#FDFFB6",
  "BioMe"              = "#BEB8D8",
  "CCPM"               = "#FA8375",
  "Estionian Biobank"  = "#81AED0",
  "Genes & Health"     = "#FBB569",
  "Genomics England"   = "#AFDD70",
  "MGBB"               = "#FCCDE4",
  "PMBB"               = "#D8D8D8",
  "UK Biobank"         = "#BC7FBA",
  "ALL"                = "#7E7E7E"
)

make_upset_plot <- function(df, xvar,
                           color_map = NULL,
                           title = NULL) {

  df <- df %>%
    group_by(id) %>%
    summarise(set = list(unique(!!sym(xvar))), .groups = "drop")

  df <- df %>%
    mutate(set_str = vapply(set, function(x) paste(sort(x), collapse = "_"), character(1)))

  # --- Helper: convert to hex ---
  to_hex <- function(cols) {
    rgb_mat <- grDevices::col2rgb(cols)
    grDevices::rgb(rgb_mat[1, ], rgb_mat[2, ], rgb_mat[3, ],
                   maxColorValue = 255)
  }

  mix_set_color <- function(set_name, base_colors) {
    parts <- unlist(strsplit(as.character(set_name), "_"))
    cols <- base_colors[parts]
    cols <- cols[!is.na(cols)]

    if (length(cols) == 0) return("grey70")

    cols <- to_hex(cols)

    if (length(cols) == 1) return(unname(cols[1]))

    mixed <- colorspace::hex2RGB(cols[1])
    for (i in 2:length(cols)) {
      mixed <- colorspace::mixcolor(0.5, mixed, colorspace::hex2RGB(cols[i]))
    }

    colorspace::hex(mixed)
  }

  sets <- unique(df$set_str)

  fill_cols <- setNames(
    vapply(sets, mix_set_color, character(1), base_colors = color_map),
    sets
  )

  p <- ggplot(df, aes(x = set, fill = set_str)) +
    geom_bar() +
    geom_text(
      stat = "count",
      aes(label = after_stat(count)),
      vjust = -0.3,
      size = 2.5
    ) +
    scale_fill_manual(values = fill_cols) +
    scale_x_upset(order_by = "freq") + 
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    # coord_cartesian(clip = "off") +
    theme_minimal() +
    labs(x = NULL, y = "Count", title = title) +
    guides(fill = "none")

  return(p)
}

create_all_plots <- function(dt_meta, file_out_append="_biobank", width_p=6, height_p=3.5)
{
    dt_meta_cauchy <- dt_meta %>% filter(source == "Cauchy")
    meta_biobanks <- dt_meta_cauchy %>%
      distinct(Region, phenotype, biobank) %>%
      mutate(id = paste(Region, phenotype, sep = "_"))

    # UpSet
    p_all <- make_upset_plot(meta_biobanks, "biobank", color_map = biobank_colors)
    pdf(paste0("Figures/upset_plot", file_out_append, ".pdf"),
      width = width_p, height = height_p)
    print(p_all)
    dev.off()
}

create_all_plots(cauchy, width_p = 13)

cauchy <- fread("../cauchy_combined_results.tsv.gz") %>% filter(Cauchy_Pvalue < (0.05/20000)) %>%
  filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
  rename(biobank = meta)
cauchy$source <- "Cauchy"

cauchy <- cauchy %>% mutate(biobank = ifelse(biobank == "ALL", "ALL", unlist(renaming_plot_biobank_list[biobank])))
create_all_plots(cauchy, file_out_append="_biobank_only", width_p = 9)

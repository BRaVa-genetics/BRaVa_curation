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
pop_colors[["non-EUR"]] <- "#1F3B73"
pop_colors[["ALL"]] <- "grey50"
pop_colors[["UKB and AoU"]] <- "#155f6f"
colors[["pLoF;damaging_missense_or_protein_altering"]] <- "#74099E"
colors[["pLoF;DM/PA"]] <- "#74099E"
colors[["DM/PA"]] <- "#FF6103"

colors_class <- c(Burden = "#ff7f00", SKAT = "#1f78b4", `SKAT-O` = "#33a02c")

meta <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_all_meta_subsets_101025.tsv.gz") %>%
# meta <- fread("../significant_assocs_from_all_meta_subsets_101025.tsv.gz") %>%
  filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
  mutate(ancestry = case_when(
    ancestry == "all" ~ "Meta",
    ancestry == "non_EUR" ~ "non-EUR",
    ancestry == "just_uk-biobank_and_all-of-us" ~ "UKB and AoU",
    TRUE ~ ancestry
  )) %>% 
  mutate(Group = case_when(
  	Group == "damaging_missense_or_protein_altering" ~ "DM/PA",
  	Group == "pLoF;damaging_missense_or_protein_altering" ~ "pLoF;DM/PA",
  	TRUE ~ Group))

cauchy <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_cauchy_assocs_from_all_meta_subsets_101025.tsv.gz") %>% 
# cauchy <- fread("../significant_cauchy_assocs_from_all_meta_subsets_101025.tsv.gz") %>%
  filter(!(Region %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))) %>%
  mutate(ancestry = case_when(
    ancestry == "all" ~ "Meta",
    ancestry == "non_EUR" ~ "non-EUR",
    ancestry == "just_uk-biobank_and_all-of-us" ~ "UKB and AoU",
    TRUE ~ ancestry
  ))

meta$source <- "Bonferroni"
cauchy$source <- "Cauchy"

combined <- bind_rows(meta, cauchy)
combined <- combined %>% mutate(ancestry = ifelse(ancestry == "Meta", "ALL", ancestry))

get_euler_weights <- function(named_list) {
  all_ids <- unique(unlist(named_list))
  sets <- names(named_list)
  binary_mat <- sapply(sets, function(set) all_ids %in% named_list[[set]])
  rownames(binary_mat) <- all_ids
  membership <- apply(binary_mat, 1, function(x) paste(sets[which(x)], collapse = "&"))
  weights <- table(membership)
  setNames(as.numeric(weights), names(weights))
}

make_upset_plot <- function(df, xvar,
                           color_map = NULL,
                           title = NULL) {

  df <- df %>%
    group_by(id) %>%
    summarise(set = list(unique(!!sym(xvar))), .groups = "drop")

  # --- Create string version ---
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
    theme_minimal() +
    labs(x = NULL, y = "Count", title = title) +
    guides(fill = "none")

  return(p)
}

make_euler_plot <- function(weights, color_map = NULL, title = NULL) {
  fit <- euler(weights, shape = "ellipse")
  fills <- NULL
  edge_cols <- TRUE
  if (!is.null(color_map)) {
    vars <- unique(unlist(strsplit(names(weights), "&")))
    fills <- list(fill = color_map[vars], alpha = 0.6)
    edge_cols <- list(col = color_map[vars], lwd = 1.5)
  }

  plot(fit,
       fills = fills,
       labels = list(font = 2),
       edges = edge_cols,
       quantities = TRUE,
       main = title)
}

create_all_plots <- function(dt_meta, file_out_append="") {
  
  if ("source" %in% names(dt_meta)) {
    dt_meta_bonf <- dt_meta %>% filter(source == "Bonferroni")
    dt_meta_cauchy <- dt_meta %>% filter(source == "Cauchy")
  } else {
    dt_meta_bonf <- dt_meta
    dt_meta_cauchy <- dt_meta
  }

  # Output Euler and UpSet: Class
  meta_class <- dt_meta_bonf %>%
    filter(ancestry == "ALL") %>%
    distinct(Region, phenotype, class) %>%
    mutate(id = paste(Region, phenotype, sep = "_"))

  venn_data_class <- split(meta_class$id, meta_class$class)
  weights_class <- get_euler_weights(venn_data_class)

  pdf(paste0("Figures/upset_plot_class", file_out_append, ".pdf"),
    width = 2.5, height = 3.5)
  print(make_upset_plot(meta_class, "class", color_map = colors_class))
  dev.off()

  pdf(paste0("Figures/euler_plot_class", file_out_append, ".pdf"),
    width = 1.5, height = 1.5)
  print(make_euler_plot(weights_class, color_map = colors_class))
  dev.off()

  # Output Euler and UpSet: Variant Group
  meta_group <- dt_meta_bonf %>%
    filter(ancestry == "ALL") %>%
    distinct(Region, phenotype, Group) %>%
    mutate(id = paste(Region, phenotype, sep = "_"))

  venn_data_group <- split(meta_group$id, meta_group$Group)
  weights_group <- get_euler_weights(venn_data_group)

  pdf(paste0("Figures/upset_plot_group", file_out_append, ".pdf"),
    width = 2.5, height = 3.5)
  print(make_upset_plot(meta_group, "Group", color_map = colors))
  dev.off()

  pdf(paste0("Figures/euler_plot_group", file_out_append, ".pdf"),
    width = 1.5, height = 1.5)
  print(make_euler_plot(weights_group, color_map = colors))
  dev.off()

  # Ancestry Subsets (ALL/EUR/non-EUR and UKB/AoU)
  meta_venn <- dt_meta_cauchy %>%
    filter(ancestry %in% c("ALL", "EUR", "non-EUR", "UKB and AoU")) %>%
    distinct(Region, phenotype, ancestry) %>%
    mutate(id = paste(Region, phenotype, sep = "_"))

  list_pops <- list(
    superpops = c("ALL", "EUR", "non-EUR"),
    ukb_aou = c("ALL", "UKB and AoU")
  )

  widths <- c(2.5, 2)
  i <- 1

  for (subset_name in names(list_pops)) {
    pops <- list_pops[[subset_name]]
    subset_df <- meta_venn %>% filter(ancestry %in% pops)
    venn_list <- split(subset_df$id, subset_df$ancestry)
    weights <- get_euler_weights(venn_list)

    # UpSet
    p <- make_upset_plot(subset_df, "ancestry", color_map = pop_colors)
    pdf(paste0("Figures/upset_plot_", subset_name, file_out_append, ".pdf"),
      width = widths[i], height = 3.5)
    print(p)
    dev.off()

    # Euler
    pdf(paste0("Figures/euler_plot_", subset_name, file_out_append, ".pdf"),
      width = 1.5, height = 1.5)
    print(make_euler_plot(weights, color_map = pop_colors))
    dev.off()

    i <- i + 1
  }

  # Final: All Ancestries
  meta_all <- dt_meta_cauchy %>%
    filter(ancestry %in% c("ALL", "AFR", "AMR", "EAS", "EUR", "SAS")) %>%
    distinct(Region, phenotype, ancestry) %>%
    mutate(id = paste(Region, phenotype, sep = "_"))

  venn_data_all <- split(meta_all$id, meta_all$ancestry)
  weights_all <- get_euler_weights(venn_data_all)

  # Euler
  pdf(paste0("Figures/euler_plot_all_ancestries", file_out_append, ".pdf"),
    width = 2.5, height = 2.5)
  print(make_euler_plot(weights_all, color_map = pop_colors))
  dev.off()

  # UpSet
  p_all <- make_upset_plot(meta_all, "ancestry", color_map = pop_colors)
  pdf(paste0("Figures/upset_plot", file_out_append, ".pdf"),
    width = 6, height = 3.5)
  print(p_all)
  dev.off()
}

create_all_plots(combined)
create_all_plots(combined %>% filter(phenotype != "Height_ALL"), file_out_append="_no_height")

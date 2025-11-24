#!/bin/Rscript
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggupset)
library(eulerr)

source("../meta_analysis_utils.r")  # must include pop_colors here
pop_colors[["non-EUR"]] <- "#1F3B73"
pop_colors[["Meta"]] <- "grey50"
pop_colors[["UKB and AoU"]] <- "#155f6f"
colors[["pLoF;damaging_missense_or_protein_altering"]] <- "#74099E"
colors[["pLoF;DM/PA"]] <- "#74099E"
colors[["DM/PA"]] <- "#FF6103"

colors_class <- c(Burden = "#ff7f00", SKAT = "#1f78b4", `SKAT-O` = "#33a02c")

meta <- fread("/well/lindgren/dpalmer/BRaVa_meta-analysis_outputs/significant_assocs_from_all_meta_subsets_101025.tsv.gz") %>%
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

get_euler_weights <- function(named_list) {
  all_ids <- unique(unlist(named_list))
  sets <- names(named_list)
  binary_mat <- sapply(sets, function(set) all_ids %in% named_list[[set]])
  rownames(binary_mat) <- all_ids
  membership <- apply(binary_mat, 1, function(x) paste(sets[which(x)], collapse = "&"))
  weights <- table(membership)
  setNames(as.numeric(weights), names(weights))
}

make_upset_plot <- function(df, xvar, color_map = NULL, title = NULL) {
  df <- df %>%
    group_by(id) %>%
    summarise(set = list(unique(!!sym(xvar))), .groups = "drop") %>%
    mutate(set_str = sapply(set, function(x) paste(sort(x), collapse = "_")))

  p <- ggplot(df, aes(x = set, fill = set_str)) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 2.5) +
    scale_x_upset() +
    theme_minimal() +
    labs(x = NULL, y = "Count", fill = xvar, title = title) +
    theme(legend.position = "none")

  if (!is.null(color_map)) {
    p <- p + scale_fill_manual(values = color_map)
  }

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
  # Output Euler and UpSet: Class
  meta_class <- dt_meta %>%
    filter(ancestry == "Meta") %>%
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
  meta_group <- dt_meta %>%
    filter(ancestry == "Meta") %>%
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

  # Ancestry Subsets (Meta/EUR/non-EUR and UKB/AoU)
  meta_venn <- dt_meta %>%
    filter(ancestry %in% c("Meta", "EUR", "non-EUR", "UKB and AoU")) %>%
    distinct(Region, phenotype, ancestry) %>%
    mutate(id = paste(Region, phenotype, sep = "_"))

  list_pops <- list(
    superpops = c("Meta", "EUR", "non-EUR"),
    ukb_aou = c("Meta", "UKB and AoU")
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
  meta_all <- dt_meta %>%
    filter(ancestry %in% c("Meta", "AFR", "AMR", "EAS", "EUR", "SAS")) %>%
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

create_all_plots(meta)
create_all_plots(meta %>% filter(phenotype != "Height_ALL"), file_out_append="_no_height")

#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(patchwork)
devtools::install_github("mkanai/rgsutil")
library(rgsutil)
library(gghighlight)
library(ggtern)
source("../meta_analysis_utils.r")

options(ggrastr.default.dpi = 300)

plt_theme <- theme_classic(base_size = 8) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "cm"),
    plot.tag = element_text(face = 'bold'),
    plot.tag.position = c(0, 1),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(size = 0)
  )


# cf: https://github.com/macarthur-lab/gnomad_lof/blob/master/R/constants.R
get_pop_colors <- function(pop = NULL) {
  color_afr <- "#941494"
  color_amr <- "#ED1E24"
  color_asj <- "coral"
  color_eas <- "#108C44"
  color_eur <- color_nfe <- "#6AA5CD"
  color_fin <- "#002F6C"
  color_mde <- "#33CC33"
  color_mid <- "#EEA9B8"
  color_oth <- "#ABB9B9"
  color_sas <- "#FF9912"
  color_oce <- "#A65628"

  color_cases <- "darkorange"
  color_controls <- "darkblue"

  pop_colors <- c(
    AFR = color_afr,
    AMR = color_amr,
    ASJ = color_asj,
    CSA = color_sas,
    EAS = color_eas,
    EUR = color_eur,
    FIN = color_fin,
    MDE = color_mde,
    MID = color_mid,
    NFE = color_nfe,
    SAS = color_sas,
    OCE = color_oce,
    cases = color_cases,
    controls = color_controls,
    Reference = "gray60",
    Remaining = "gray60"
  )

  if (!is.null(pop)) {
    return(pop_colors[pop])
  }
  return(pop_colors)
}


# cf. https://github.com/atgu/ukbb_pan_ancestry/blob/master/plot_ukbb_pca.R
plot_pca <- function(dataset, first_pc, second_pc, color_pop,
  xlim = NULL, ylim = NULL, images = NULL) {
    pc_biplot <-
      dplyr::arrange(dataset, !!as.symbol(color_pop)) %>%
      ggplot(aes_string(x = first_pc, y = second_pc, color = color_pop)) +
      images +
      ggrastr::rasterize(geom_point(alpha = 0.2)) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      plt_theme +
      scale_color_manual(values = get_pop_colors(),
                         name = "Population",
                         na.value = "grey50") +
      coord_cartesian(xlim = xlim, ylim = ylim)
    return(pc_biplot)
  }


plot_pca_density <- function(dataset, first_pc, second_pc,
  xlim = NULL, ylim = NULL) {
    pc_biplot <-
      ggplot(dataset, aes_string(x = first_pc, y = second_pc)) +
      geom_hex(bins = 50) +
      plt_theme +
      scale_fill_gradientn(
        trans = "log",
        name = "Count",
        colours = rev(RColorBrewer::brewer.pal(5, "Spectral"))
      ) +
      coord_cartesian(xlim = xlim, ylim = ylim)
    return(pc_biplot)
  }


annotate_image <- function(path, greyscale = FALSE, offset_x = 22,
  offset_y = 12) {
  image = magick::image_read(path) %>%
    magick::image_crop(
      sprintf(
        "%dx%d+%d+%d",
        1800 - offset_x,
        1800 - offset_x - offset_y,
        offset_x,
        offset_y
      )
    )
  if (greyscale) {
    image = magick::image_convert(image, colorspace = "gray")
  }
  return(list(annotation_raster(as.raster(image),-Inf, Inf,-Inf, Inf)))
}


or_missing <- function(predicate, value) {
  if (predicate) {
    value
  } else {
    NULL
  }
}


metadata <- rgsutil::read_gsfile("gs://covid19-hg-public/pca_projection/gnomad_meta_hgdp_tgp_v1.txt")
reference <- rgsutil::read_gsfile("gs://gbmi-public/hgdp_tgp_pca_gbmi_snps_scores.txt.bgz") %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(pop = hgdp_tgp_meta.Genetic.region, study = project_meta.title) %>%
  dplyr::mutate(pop = ifelse(pop == "CSA", "SAS", pop))

# This defines the plotting region
plot_pcs <- paste0("PC", seq(10))
reference_range <-
  purrr::map(plot_pcs, function(pc) {
    range(reference[[pc]])
  }) %>%
  magrittr::set_names(plot_pcs)

cohort_files <- list(
  `all-of-us` = "gs://brava-meta-upload-all-of-us/pca/all-of-us.Lu.projected.pca.tsv.gz",
  bbj = "gs://brava-meta-upload-bbj/pca_projection/bbj.Sonehara.20250215.projected.pca.tsv.gz",
  biome = "gs://brava-meta-upload-biome/projected_pc_plots_20250417/BioMe.projected.pca.tsv.gz",
  ccpm = "gs://brava-meta-upload-ccpm/brava-results/pca_projections/CCPM.projected.pca.tsv.gz",
  egcut = "gs://brava-meta-upload-egcut/EstBB.projected.pca.tsv.gz",
  `genes-and-health` = "gs://brava-meta-upload-genes-and-health/pca_projections/GNH.GKalantzis.20250211.projected.pca.tsv.gz",
  pmbb = "gs://brava-meta-upload-pmbb/pca_projections/PMBB.projected.pca.tsv.gz",
  mgbb = "gs://brava-meta-upload-mgbb/20250304/MGBB.KOYAMA.20250304.projected.pca.tsv.gz",
  `uk-biobank` = "gs://brava-meta-upload-uk-biobank/uk-biobank.palmer.projected.pca.tsv.gz"
)

read_gsfile_and_name <- function(file) {
  dt <- rgsutil::read_gsfile(file)
  s <- gsub("^.*/([A-Za-z0-9-]+)\\..*", "\\1", file)
  for (d in names(file_check_information$dataset)) {
    if (s %in% file_check_information$dataset[[d]]) {
      s <- d
    }
  }
  s <- renaming_plot_biobank_list[[s]]
  dt <- dt %>% mutate(study = s)
  return(dt %>% filter(pop != "Reference"))
}

cohort <- purrr::map_dfr(cohort_files, read_gsfile_and_name) %>%
  # plots only
  dplyr::bind_rows(
    tibble::tibble(study = renaming_plot_biobank_list[["gel"]], pop = "dummy", PC1 = reference$PC1[1], PC2 = reference$PC2[1])
  )

# Check counts make sense
table(cohort$study)

combined <- dplyr::bind_rows(
  dplyr::select(reference, dplyr::starts_with("PC"), pop, study) %>%
    dplyr::mutate(study = "Reference"),
  dplyr::select(cohort, dplyr::starts_with("PC"), pop, study)) %>%
  dplyr::mutate(
    study = factor(study, levels = c("Reference", sort(unique(cohort$study))))
  )

reference_plt <- matrix(paste0("PC", seq(10)), ncol = 2, byrow = TRUE) %>%
  magrittr::set_colnames(c("x", "y")) %>% tibble::as_tibble() %>%
  purrr::pmap(function(x, y) {
    plot_pca(reference, x, y, color_pop = "pop")
  }) %>% append(list(patchwork::guide_area())) %>%
  purrr::reduce(`+`) + patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = "a")

cowplot::save_plot(
  "Figures/reference_plot.pdf",
  reference_plt,
  base_height = 4.5,
  base_width = 7.2,
  device = cairo_pdf
)

studies <- sort(unique(combined$study))
combined <- combined %>% filter(pop != "Remaining")
ncol <- 5
padding_percent_plot <- 20
scaling_range_plot <- 1 + padding_percent_plot/100

# Plot pairs of PCs, including the image to plot over the top
for (pc in c(1,3,5,7,9))
{
  first_pc <- paste0("PC", pc)
  second_pc <- paste0("PC", pc+1)
  cat("Creating biplots for PCs...", first_pc, "versus", second_pc, "\n")
  plt <- purrr::map(
    seq_along(studies),
    function(i) {
      study <- studies[i]
      data <- dplyr::filter(combined, study  == .env$study) %>% 
        dplyr::arrange(pop)
      reference <- dplyr::filter(combined, study  == "Reference") %>%
        dplyr::arrange(pop)
      color_pop <- "pop"
      xlim <- reference_range[[first_pc]] * scaling_range_plot
      ylim <- reference_range[[second_pc]] * scaling_range_plot
      is_first_row <- (i %% ncol) == 1
      is_last_col <- ((i - 1) %/% ncol) == (length(studies) %/% ncol)

      if (study == "Genomics England") {
        images <- annotate_image(paste0(
          "gel_kousathanas_BRaVa_pca_projection_plots_20250217/",
          "gel.pca_projection.projected.pca.ancestry.no.ref.PC",
          pc, "-", pc+1, ".png"))
      } else {
        images <- NULL
      }

      ggplot() + ggrastr::rasterize(
        list(or_missing(
          study != "Reference",
          geom_point(aes_string(x = first_pc, y = second_pc), alpha = 0.2,
            size = 0.1, color = "grey50", data = reference)),
          geom_point(aes_string(x = first_pc, y = second_pc, color = color_pop),
            alpha = 0.2, size = 0.1, data = data)
          )
        ) + images + 
        guides(color = guide_legend(ncol = 2, override.aes = list(alpha = 1))) +
        plt_theme +
        theme(
          plot.title = element_text(size = 6),
          legend.title = element_text(size = 8),
          legend.key.width = unit(2, "mm"), legend.key.height = unit(5, "mm")) +
        scale_color_manual(
          values = get_pop_colors(),name = "Population", na.value = "grey50") +
        scale_x_continuous(
          labels = scales::number_format(accuracy = 0.01, drop0trailing=TRUE)) +
        scale_y_continuous(
          labels = scales::number_format(accuracy = 0.01, drop0trailing=TRUE)) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        or_missing(
          study != "Reference", theme(legend.position = "none")) +
        or_missing(
          !is_last_col, theme(
            axis.title.x = element_blank(), axis.text.x = element_blank())) +
        or_missing(
          !is_first_row, theme(
            axis.title.y = element_blank(), axis.text.y = element_blank())) +
        labs(title = study)
    }

  ) %>%
    append(list(patchwork::guide_area())) %>%
    purrr::reduce(`+`) +
    patchwork::plot_layout(ncol = ncol, guides = 'collect')

  cowplot::save_plot(paste0("Figures/PCA_", first_pc, "_", second_pc,
    "_per_cohort.pdf"), plt, base_height = 5.6, base_width = 7.2,
  device = cairo_pdf)
}

# Ternery plots
dat.tern <- dplyr::select(combined, PC1, PC2, PC3)  %>%
  dplyr::mutate(
    PC1 = -PC1,
    PC2 = -PC2,
    PC3 = -PC3,
    PC1 = PC1 + (min(PC1, na.rm = TRUE) * -1),
    PC2 = PC2 + (min(PC2, na.rm = TRUE) * -1),
    PC3 = PC3 + (min(PC3, na.rm = TRUE) * -1)
  ) %>% as.matrix() %>% prop.table(., 1)  %>%
  as.data.frame() %>% dplyr::bind_cols(dplyr::select(combined, study, pop))

plt <- dplyr::mutate(dat.tern,
    panel = factor(
      ifelse(study == "Reference", "Reference", "BRaVa cohorts"),
    levels = c("Reference", "BRaVa cohorts"))) %>% 
  dplyr::filter(study != "Reference") %>%
  ggplot(aes(x = PC1, y = PC2, z = PC3, colour = pop)) +
  ggtern::coord_tern() + 
  ggrastr::rasterize(geom_point(size = 0.5, alpha = 0.1)) +
  guides(color = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
  scale_color_manual(values = get_pop_colors(), name = "Population",
    na.value = "grey50") + plt_theme +
  theme(legend.position = "bottom", 
        legend.justification = "bottom",
        legend.key.width = unit(0.3, "cm")) +
  ggtern::theme_showarrows() +  ggtern::theme_noticks() +
  ggtern::theme_nolabels() + ggtern::theme_notitles() +
  theme(strip.text = element_text(size = 12))

ggtern::ggsave("Figures/ternery.pdf", plt, height = 3.6, width = 3.6)
system("mv Figures/PCA_PC1_PC2_per_cohort.pdf Figures_main_text/Figure2.pdf")
system("mv Figures/PCA_PC3_PC4_per_cohort.pdf Figures_supplement/FigureS1.pdf")
system("mv Figures/PCA_PC5_PC6_per_cohort.pdf Figures_supplement/FigureS2.pdf")

library(googlesheets4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

sheet_url <- "https://docs.google.com/spreadsheets/d/1LMtyhOIWgyPfs_UGv2CA6yh8jnXV7gCjN_5zRSi7ORw/edit?gid=179641139#gid=179641139"

pdf(width=4, height=4, file="confusion.pdf")
for (s in c("Table S11", "Table S12", "Table S13", "Table S14")) {
  dat <- read_sheet(sheet_url, sheet = s)

  lit_col  <- "Literature verdict"
  open_col <- "Open targets verdict"

  lvl <- c("Not found", "Hypothesized", "Existing", "Established")

  dat2 <- dat %>%
    transmute(
      literature = factor(.data[[lit_col]], levels = lvl),
      open_targets = factor(.data[[open_col]], levels = lvl)
    )

  plot_df <- dat2 %>%
    count(literature, open_targets, name = "n") %>%
    complete(literature = lvl, open_targets = lvl, fill = list(n = 0)) %>%
    mutate(
      literature = factor(literature, levels = lvl),
      open_targets = factor(open_targets, levels = lvl),
      text_color = ifelse(n >= max(n) * 0.45, "white", "#1a1a1a")
    )

  p <- ggplot(plot_df, aes(x = open_targets, y = literature, fill = n)) +
    geom_tile(color = "white", linewidth = 1) +
    geom_text(aes(label = n, color = text_color), size = 5, fontface = "bold", show.legend = FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(
      colors = c(
        "#f7fbff",
        "#deebf7",
        "#c6dbef",
        "#9ecae1",
        "#6baed6",
        "#3182bd",
        "#08519c"
      )
    ) +
    coord_equal() +
    labs(
      title = "Confusion Matrix (Lit vs OT)",
      x = "Open Targets",
      y = "Literature",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      legend.position = "none"
    )

  print(p)
}

dev.off()
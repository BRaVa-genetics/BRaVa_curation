library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)

# Results to exclude
# here on the cluster /well/lindgren/dpalmer/BRaVa_meta-analysis_inputs/inflation_summaries.tsv.gz
dt_inflation <- fread("../inflation_summaries.tsv.gz")
dt_inflation <- dt_inflation %>% filter(Group == "synonymous") %>% filter(max_MAF != 0.01, lambda_value < 1.3)

pdf("Figures/inflation.pdf", width=10, height=10)

# Ensure variables are factors
dt <- dt_inflation
dt[, lambda_type := factor(lambda_type)]
dt[, max_MAF := factor(max_MAF, levels = sort(unique(max_MAF)))]

# Convert lambda_type to label-parsed format
dt[, lambda_type := factor(lambda_type)]
dt[, lambda_type := factor(sapply(as.character(lambda_type), function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) == 3) {
    sprintf("lambda[%s]^'%s'", parts[2], parts[3])
  } else if (length(parts) == 2) {
    sprintf("lambda[%s]", parts[2])
  } else {
    "lambda"
  }
}))]

# Convert max_MAF to label-parsed format (e.g., "10^{-4}")
dt[, max_MAF := factor(as.character(max_MAF))]
dt[, max_MAF := factor(sapply(as.character(max_MAF), function(x) {
  if (grepl("^1e-", x)) {
    sprintf("10^{-%s}", sub("1e-", "", x))
  } else {
    x
  }
}))]

# Plot for each ancestry
for (anc in unique(dt$ancestry)) {
  dt_anc <- dt[ancestry == anc]

  p <- ggplot(dt_anc, aes(x = dataset, y = lambda_value)) +
    geom_boxplot() +
    facet_grid(
      rows = vars(max_MAF),
      cols = vars(lambda_type),
      labeller = label_parsed
    ) +
    theme_bw() +
    labs(
      title = paste("Ancestry:", anc),
      x = "Dataset",
      y = expression(lambda~"Value")
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10)
    )

  print(p)
}
library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

source("../meta_analysis_utils.r")

# Create the plots of prevalence and count, by-disease
sampling_colors <- c(
  "Hospital/Health center–based" = "#2b2b2b",
  "Mixed" = "#999999",
  "Population-based" = "#e0e0e0"
)

make_trait_panel <- function(df, trait_label) 
{
	df <- df %>% mutate(Biobank = factor(Biobank, levels = unique(Biobank)))
	d <- df %>% filter(phenotype == trait_label)    
	trait_name <- renaming_phenotype_list[trait_label]

	d_prev <- d %>% group_by(Biobank) %>%
		summarise(
			Prevalence = sum(as.integer(n_cases)) /
			(sum(as.integer(n_cases)) + sum(as.integer(n_controls))),
		sampling_strategy = first(sampling_strategy)) %>%
		arrange(Prevalence) %>%
		mutate(Biobank = factor(Biobank, levels = Biobank)) 
	levels_order <- levels(d_prev$Biobank)

	# Apply to d2
	d[, Biobank := factor(Biobank, levels = levels_order)]

	p1 <- ggplot(d_prev,
		aes(x = Prevalence*100, y = Biobank, fill = sampling_strategy)) +
	geom_col() +
	scale_fill_manual(values = sampling_colors) +
	labs(x = "Prevalence (%)", y = NULL) +
	theme_minimal() +
	theme(legend.position = "none",
	      plot.title = element_text(hjust = 0.5))

	p2 <- ggplot(d, aes(x = n_cases, y = Biobank, fill = ancestry)) +
	geom_col() +
	scale_fill_manual(values = pop_colors) +
	labs(x = "Number of cases", y = NULL) +
	theme_minimal() +
	theme(legend.position = "none")

  # Combine two panels horizontally
  (p1 + p2) + 
  plot_annotation(title = trait_name)  # or use any string
}

dt_info <- fread("../manuscript_tables/Tables/table_for_creation_of_prevalence_plots.tsv") %>%
	filter(phenotype %in% 
		c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
		"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
		"CRF", "CRP", "EFRMB", "FemInf","Gout", "HDLC", "Height", "HF",
		"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
		"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
		"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
		"WHRBMI", "AlcCons", "ALT", "HipRep"))

# Loop over traits
all_traits <- unique(dt_info$phenotype)
plot_list <- lapply(all_traits, function(tr) make_trait_panel(dt_info, tr))

pdf("Figures/trait_prevalences.pdf", width=7, height=2)
lapply(plot_list, print)
dev.off()

# Create pngs as well to paste into the google doc
for (i in 1:length(all_traits)) {
	ggsave(plot_list[[i]],
		filename=paste0("Figures/", all_traits[i], "_prevalence.png"),
		width=7, height=2)
}


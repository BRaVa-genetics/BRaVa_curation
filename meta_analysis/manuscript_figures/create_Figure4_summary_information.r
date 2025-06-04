library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(latex2exp)

# Loop over the data files that contain the case and control counts
files <- dir(path="../../data/meta_analysis/gcloud", full.names=TRUE, recursive=TRUE)
files <- gsub("^.*/", "", grep("cleaned", files, value=TRUE))

source("../meta_analysis_utils.r")

file_info <- rbindlist(lapply(files, extract_file_info), fill=TRUE)
file_info <- file_info %>% filter(phenotype %in% 
	c("AFib", "AMD", "AST", "Asth", "BenCervUterNeo", "BenIntNeo",
	"BMI", "BreastCanc", "CAD", "CervCanc", "ColonRectCanc", "COPD",
	"CRF", "EFRMB", "FemInf","Gout", "HDLC", "HF",
	"HTN", "IBD", "IFHern", "ILDSarc", "LDLC", "MatHem", "NonRheuValv",
	"PAD", "Pancreat", "PeptUlcer", "Psori", "RheumArth", "RheumHeaDis",
	"Stroke", "T2Diab", "TChol", "TG", "Urolith", "VaricVeins", "VTE",
	"WHRBMI", "AlcCons", "ALT", "HipRep"))

renaming_phenotype_list[unique(file_info$phenotype)]

ancestry_information <- file_info %>% group_by(dataset, ancestry) %>%
	summarise(n=max(as.integer(n), (as.integer(n_cases) + as.integer(n_controls)), na.rm=TRUE), .groups='drop') %>%
	complete(dataset, ancestry, fill = list(n=0))

ancestry_information <- ancestry_information %>% group_by(dataset) %>% mutate(proportion = n/sum(n))

# Rename
ancestry_information <- ancestry_information %>%
	mutate(dataset = unlist(renaming_plot_biobank_list[dataset]))

# Plot the results

# First the counts
pdf(width = 3.7, height = 3, file = "Figures/ancestry_summary.pdf")
p <- ggplot(ancestry_information, aes(x=dataset, y=n, fill=ancestry)) + 
	geom_bar(stat= "identity", position = "dodge") + theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  	scale_y_continuous(labels = scales::comma) +
  	scale_fill_manual(values=pop_colors) + 
  	labs(x = NULL, y = "Number of samples", fill = "Genetic\nancestry")
print(p)

p <- ggplot(ancestry_information, aes(x=dataset, y=n, fill=ancestry)) + 
	geom_bar(stat= "identity", position = "stack") + theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  	scale_y_continuous(labels = scales::comma) +
  	scale_fill_manual(values=pop_colors) + 
  	labs(x = NULL, y = "Number of samples", fill = "Genetic\nancestry")
print(p)

# Second the proportions
p <- ggplot(ancestry_information, aes(x=dataset, y=n, fill=ancestry)) + 
	geom_bar(stat= "identity", position = "fill") + theme_minimal() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  	scale_y_continuous(labels = scales::comma) +
  	scale_fill_manual(values=pop_colors) + 
  	labs(x = NULL, y = "Proportion of samples", fill = "Genetic\nancestry")
print(p)

# Next, let's determine how many phenotypes are represented for each of the cohorts in the meta-analysis
n_phenotypes_summary <- file_info %>% filter(type == "gene", (is.na(n) & ((n_cases > 100) | (n_controls > 100))) | !is.na(n)) %>%
# n_phenotypes_summary <- unique(n_phenotypes_summary %>% select(dataset, phenotype, ancestry)) %>%
	group_by(dataset, ancestry, phenotype) %>% slice(1) %>% group_by(dataset, ancestry) %>% summarise(count = n(), .groups = 'drop') %>%
	complete(dataset, ancestry, fill = list(count=0))

p <- ggplot(n_phenotypes_summary %>% 
	mutate(
		dataset = unlist(renaming_plot_biobank_list[dataset])),
		aes(x=dataset, y=count, fill=ancestry)) + 
	theme_minimal() + geom_bar(stat="identity", position="dodge", width=0.7) + 
	scale_fill_manual(values=pop_colors) + 
	labs(x = NULL, y = "Number of phenotypes", fill = "Genetic\nancestry") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()
 
n_phenotypes_summary <- file_info %>% filter(type == "gene", (is.na(n) & ((n_cases > 100) | (n_controls > 100))) | !is.na(n)) %>% 
	group_by(dataset, phenotype) %>% slice(1) %>% group_by(dataset) %>% summarise(count = n(), .groups = 'drop')

pdf(width = 2.5, height = 3, file = "Figures/phenotype_counts.pdf")

p <- ggplot(n_phenotypes_summary %>% 
	mutate(
		dataset = unlist(renaming_plot_biobank_list[dataset])),
		aes(x=dataset, y=count, fill=dataset)) + 
	theme_minimal() + geom_bar(color = "darkgrey", stat="identity", position="dodge", width=0.7, show.legend = FALSE) + 
	labs(x = NULL, y = "Number of phenotypes") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_fill_brewer(palette = "Set3")

print(p)

dev.off()

run <- ancestry_information %>% group_by(dataset) %>% summarise(total_n = sum(n)) %>% mutate(run = TRUE)
# Add in estimates of number of samples in the remaining cohorts (at present)
additional_pledged_samples <- data.table(
	dataset = c("ALSPAC", "CKB", "CCPM", "deCODE", "ESTBB", "DanRaV", "Qatar", "VIKING genes"),
	total_n = c(3000, 10000, 40000, 100000, 3000, 8000, 10000, 4000)
	) %>% mutate(run = FALSE)
dt_pledged <- rbind(run, additional_pledged_samples) %>% group_by(run) %>% summarise(total_n = sum(total_n))

pdf(width = 4.5, height = 4, file = "Figures/BRaVa_pledged.pdf")

p <- ggplot(dt_pledged, aes(y = total_n, x = factor(1), fill = as.character(run))) + 
  theme_minimal() + 
  geom_bar(stat = "identity", position = "stack", width = 0.5, show.legend = FALSE) + 
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = "Number of samples") +
  scale_fill_manual(values = c("TRUE" = "#E4B002", "FALSE" = "#024CFF")) +
  scale_x_discrete(labels = "BRaVa") +  # Add custom x-axis label
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10 * 1.5)  # Make x-axis label twice as big
  )

print(p)

dev.off()


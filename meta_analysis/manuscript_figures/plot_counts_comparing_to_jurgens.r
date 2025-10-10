library(googlesheets4)
library(dplyr)
library(ggplot2)

extract_phenotype_counts <- function() {
	dt <- read_sheet("https://docs.google.com/spreadsheets/d/1i3cfSFOiqHHDT6fHixe6TO7zbwbINBEj_R5r96ML7LE/edit?gid=481679696#gid=481679696",
		sheet="BRaVa_pilot_disease_endpoint_numbers", skip=1)
	names(dt) <- gsub("...[0-9]+", "", names(dt))
	names(dt)[6:10] <- paste0(names(dt)[6:10], "_jurgens")
	return(dt)
}

dt <- extract_phenotype_counts()
dt <- dt %>% filter(!is.na(`N cases_jurgens`))
names(dt)[6:10] <- gsub("_jurgens", "", names(dt)[6:10])
dt[, 6] <- dt[, 1]
dt <- rbind(dt[,1:5] %>% mutate(cohort="BRaVa"), dt[, 6:10] %>% mutate(cohort="Jurgens et al."))

pdf("BRaVa_jurgens_counts.pdf", width=10, height=5)
p <- ggplot(dt, aes(x = `Binary phenotype`, y = `N cases`, fill = cohort)) +
	geom_bar(stat = "identity", position = "dodge") +
	labs(title = "Case counts",
		x = "Phenotype",
		y = "N cases") +
	theme_minimal() + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
		axis.title.x = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 50))
print(p)
p <- ggplot(dt, aes(x = `Binary phenotype`, y = `Neff`, fill = cohort)) +
	geom_bar(stat = "identity", position = "dodge") +
	labs(title = "N effective estimates",
		x = "Phenotype",
		y = "N effective") +
	theme_minimal() + 
	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
		axis.title.x = element_blank(), plot.margin = margin(5.5, 5.5, 5.5, 50))
print(p)
dev.off()
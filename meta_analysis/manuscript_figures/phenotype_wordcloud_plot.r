# Scatter of counts is generated in 'counts_of_number_of_individuals_per_phenotype.r'

# Let's get wordclouds of the analysed traits, sized by sample size
# Then another one of the pilot traits, again sized by available sample size
library(ggwordcloud)
# Use the googlesheets package to extract the data
library(googlesheets4)

source("~/Repositories/BRaVa_curation/phenotypes/BRaVa_phenotypes_utils.r")
dt_nominated <- extract_BRaVa_pilot_phenotypes_ukb_counts(pilot_only=FALSE) %>% 
	mutate(case_control = ifelse(case_control == "N/A (quantitative trait)", FALSE, TRUE))
dt_nominated <- dt_nominated %>% mutate(word = unlist(renaming_phenotype_list[dt_nominated$phenotypeID]), freq=count) %>% select(word, freq, case_control)
dt_pilot <- extract_BRaVa_pilot_phenotypes_ukb_counts() %>% 
	mutate(case_control = ifelse(case_control == "N/A (quantitative trait)", FALSE, TRUE))
dt_pilot <- dt_pilot %>% mutate(word = unlist(renaming_phenotype_list[dt_pilot$phenotypeID]), freq=count) %>% select(word, freq, case_control)

pdf(file="Figures/wordcloud.pdf", width=5, height=5)
p <- ggwordcloud2(dt_pilot %>% 
	mutate(freq = sqrt(sqrt(ifelse(case_control, freq, 0.1*freq)))) %>%
	select(word, freq), size = 0.2, rotateRatio = 0.6, grid_margin=0.1, grid_size=2)
print(p)
p <- ggwordcloud2(dt_nominated %>% 
	mutate(freq = sqrt(sqrt(ifelse(case_control, freq, 0.1*freq)))) %>%
	select(word, freq), size = 0.2, rotateRatio = 0.6, grid_margin=0.1, grid_size=2)
print(p)
dev.off()

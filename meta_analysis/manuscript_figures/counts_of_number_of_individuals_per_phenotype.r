library(data.table)
library(dplyr)

# Location of the data files
source("meta_analysis_utils.r")
source("../QC/utils/pretty_plotting.r")

extract_counts_from_files <- function(files_folder, regexp,
	greater_than_one=FALSE) {
	files <- dir(path=files_folder, recursive=TRUE)
	gene_files <- grep(regexp, files, value=TRUE)
	gene_files <- gsub(".*/", "", gene_files)
	info_list <- list()
	for (file in gene_files) {
		info_list[[file]] <- extract_file_info(file)
		setDT(info_list[[file]], check.names=TRUE)
	}
	info <- rbindlist(info_list, fill=TRUE)

	# Determine the counts per phenotype
	dt_summary_binary <- info %>% filter(binary) %>% 
		group_by(phenotype, ancestry) %>% 
		summarise(
			total_n_cases = sum(as.integer(n_cases)),
			total_n_controls = sum(as.integer(n_controls)),
			n_biobanks=length(unique(dataset))
			)# %>% mutate(total_n = total_n_cases + total_n_controls)
	
	dt_summary_cts <- info %>% filter(!binary) %>% 
		group_by(phenotype, ancestry) %>% 
		summarise(total_n = sum(as.integer(n)),
			n_biobanks=length(unique(dataset)))
	
	dt_summary <- rbind(
		data.table(dt_summary_cts),
		data.table(dt_summary_binary), fill=TRUE)

	if (greater_than_one) {
		dt_summary <- dt_summary %>% filter(n_biobanks > 1)
	}

	return(dt_summary)
}

dt_summary_pop <- extract_counts_from_files(
	"../data/meta_analysis/gcloud",
	"cleaned.*\\.gene.*", greater_than_one = FALSE)

dt_summary <- dt_summary_pop %>% group_by(phenotype) %>% 
	summarise(
			total_n_cases = sum(total_n_cases),
			total_n_controls = sum(total_n_controls),
			total_n = sum(total_n))

dt_summary_UKB <- extract_counts_from_files(
	"../data/meta_analysis/gcloud",
	"cleaned.*uk-biobank.*EUR.*\\.gene.*",
	greater_than_one=FALSE
	)

dt_summary_UKB <- dt_summary_UKB %>% rename(
	total_n_ukb = total_n,
	total_n_cases_ukb = total_n_cases ,
	total_n_controls_ukb = total_n_controls) %>% select(-n_biobanks)

dt_summary_UKB <- data.table(dt_summary_UKB, key="phenotype")
dt_summary <- data.table(dt_summary, key="phenotype")

dt_summary <- merge(dt_summary, dt_summary_UKB)
dt_summary <- dt_summary %>% mutate(
	proportion_increase_n = total_n/total_n_ukb,
	proportion_increase_n_cases = total_n_cases/total_n_cases_ukb,
	proportion_increase_n_controls = total_n_controls/total_n_controls_ukb)

dt_ukb <- melt(dt_summary %>% 
	select(-c("total_n_cases", "total_n_controls", "total_n")) %>% rename(
	`Cases` = total_n_cases_ukb,
	`Controls` = total_n_controls_ukb,
	`Continuous count` = total_n_ukb),
	id.vars = "phenotype", 
	measure.vars = c("Cases", "Controls", "Continuous count"),
	variable.name = "type",
	value.name = "UK Biobank, EUR")

# Melt the data.table for Meta columns
dt_meta <- melt(dt_summary %>% rename(
	`Cases` = total_n_cases,
	`Controls` = total_n_controls,
	`Continuous count` = total_n), id.vars = "phenotype", 
    measure.vars = c("Cases", "Controls", "Continuous count"),
    variable.name = "type",
    value.name = "meta-analysis v2")

# Add a column to distinguish between UKB and Meta
dt_ukb[, source := "UK Biobank, EUR"]
dt_meta[, source := "meta-analysis v2"]

setkeyv(dt_meta, c("phenotype", "type"))
setkeyv(dt_ukb, c("phenotype", "type"))

dt_long <- merge(dt_meta, dt_ukb) %>% filter(phenotype != "HipRep")
dt_long <- dt_long %>% mutate()
# Combine counts
p <- create_pretty_scatter(dt_long, aes(x=`UK Biobank, EUR`, y=`meta-analysis v2`, col=factor(type)),
	x_label="UK Biobank (EUR) count", y_label="BRaVa pilot count") + theme_minimal()
pdf("count_scatter.pdf", width=8, height=4)
p <- p + 
	scale_y_continuous(labels = scales::comma) +
	scale_x_continuous(labels = scales::comma) +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	geom_abline(intercept = 0 , slope = 1, colour = 'grey40', linetype = "dashed")
print(p)
dev.off()

pdf("count_scatter_square.pdf", width=4, height=4)
p <- p + labs(title = NULL, subtitle = NULL) + theme(
    legend.position = c(.05, 1.05),
    legend.justification = c("left", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
    )
print(p)
dev.off()

# Create histograms of the scale up
dt_hist <- melt(dt_summary %>% rename(
	`Multiplier of cases` = proportion_increase_n_cases,
	`Multiplier of controls` =  proportion_increase_n_controls,
	`Multiplier of total` = proportion_increase_n), id.vars = "phenotype", 
    measure.vars = c("Multiplier of cases", "Multiplier of controls", "Multiplier of total"),
    variable.name = "type",
    value.name="Factor")

p <- create_pretty_hist(dt_hist, aes(x=Factor), x_label = "Factor", binwidth=0.2, print_p=FALSE)
p <- p + facet_wrap(~type)
pdf("factor_hist.pdf", width=6, height=3)
print(p)
dev.off()

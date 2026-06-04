library(data.table)
library(dplyr)
library(googlesheets4)

dt <- read_sheet("https://docs.google.com/spreadsheets/d/1oHuzxK0pLuekCncAUWPnEYL5h-pihti-_DpUkCdkMC0/edit#gid=1616441641", sheet="gnomAD 450k, Jul 2023", skip=2)
dt <- dt[1:which(dt$chrom == "Total"),]
dt <- dt[-which(dt$chrom == "Y"),]

for (n in names(dt)) {
	dt[n] <- unlist(dt[n])
}

# Remove the columns containing results we don't require
dt <- dt[, -grep("Cost|cost|£|GiB|Time|time|VCFs", names(dt))]
names(dt)[grep("variants", names(dt))] <- c("Raw count", "For sample metrics", "Following MAD")

fwrite(dt, file = "../../site_tables/450k/summary/00_variant_count_summary.tsv", sep='\t')

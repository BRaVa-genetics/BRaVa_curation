# Sample sizes
biobanks <- c(
		"All of Us",
		"ALSPAC",
		"BioMe",
		"Biobank Japan",
		"China Kadoorie Biobank",
		"CCPM",
		"deCODE ",
		"Estonian Biobank",
		"Indiana biobank",
		"Dan-RaV",
		"Genes & Health",
		"Genomics England",
		"Mass General Brigham Biobank",
		"Mexico City Longitudinal Study",
		"Million Veterans Program",
		"Mount Sinai Million",
		"Penn Medicine Biobank",
		"PRECISE",
		"Qatar Genome",
		"Tohoku Medical Megabank",
		"UK Biobank",
		"VIKING Genes",
		"Michigan Genomics Initiative",
		"Total")
sample_sizes <- data.table(
	Biobank = rep(biobanks, 2),
	`Sample size` = c(
		c(
		250000,
		3000,
		31000,
		10000,
		10000,
		50000,
		100000,
		3000,
		0,
		8000,
		45000,
		85000,
		45000,
		150000,
		0,
		0,
		45000,
		100000,
		10000,
		150000,
		500000,
		4000,
		600,
		1599600),
		c(1000000,
		3000,
		55000,
		200000,
		500000,
		450000,
		100000,
		200000,
		150000,
		8000,
		65000,
		85000,
		145000,
		150000,
		150000,
		1000000,
		100000,
		1000000,
		100000,
		150000,
		500000,
		4000,
		80000,
		6195000)),
	Timescale = c(rep("Now", length(biobanks)), rep("~ 5 years", length(biobanks)))
)
sample_sizes <- sample_sizes %>% mutate(
	Timescale = factor(Timescale, levels = c("Now", "~ 5 years")),
	Biobank = factor(Biobank, levels = biobanks))

ggplot(sample_sizes %>% filter(Biobank != "Total"), aes(x = Timescale, y=`Sample size`, col=Biobank)) +
geom_line(aes(group=Biobank)) + theme_minimal() +
scale_y_continuous(labels = scales::comma)

pdf("sample_sizes.pdf", width=10, height=4)
ggplot(sample_sizes %>% filter(Biobank != "Total"), aes(x = Timescale, y=`Sample size`, fill=Biobank)) +
geom_bar(stat= "identity", position = "stack") + theme_minimal() +
scale_y_continuous(labels = scales::comma) 
dev.off()

ggplot(sample_sizes %>% filter(Biobank != "Total"), aes(x = Timescale, y=`Sample size`, col=Biobank)) +
geom_line(aes(group=Biobank)) + theme_minimal() +
scale_y_continuous(labels = scales::comma) + geom_line(data=sample_sizes %>% filter(Biobank == "Total"), aes(x=Timescale, y=`Sample size`, group=Biobank), col='Black')


for (i in 1:length(file_paths)) {
	meta_list <- fread(
		paste0("~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/",
			file_root[i], "_figure_4.tsv.gz"))
	meta_list <- meta_list %>% mutate(phenotype_category = unlist(phenotype_broad_categories[phenotype]))
	# Create the Burden, SKAT, and SKAT-O versions
	# Do the same thing, but split by case-control vs cts (way more power for cts).
	for (cc in c(TRUE, FALSE)) {
		cat(ifelse(cc, "case control\n", "continuous\n"))
		meta_list_tmp <- meta_list %>% filter(case_control == cc)
		p <- make_manhattan_plot(meta_list_tmp$chromosome_name,
			meta_list_tmp$start_position,
			meta_list_tmp$Pvalue,
			threshold=1000, significance_T = 6.7e-7,
			label=meta_list_tmp$external_gene_name, 
			colour_1 = "#6583E6",
			colour_2 = "#384980")
		threshold <- ifelse(cc, 10, 50)
		p$p <- p$p + geom_label_repel(
			data = unique(subset(p$dt, y > threshold) %>% group_by(labels) %>% 
				filter(y == max(y))) %>% ungroup(),
			size = 2, aes(label=labels),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50')
		width <- 230
		height <- 100
		scaling <- 1
		file <- paste0(file_root[i], "_",
			ifelse(cc, "unique_case_control", "unique_cts"))
		ggsave(paste0(file, '.jpg'), p$p, width=width*scaling,
			height=height*scaling, units='mm')
		width <- 150
		make_gene_manhattan_category_plot(meta_list_tmp, buffer=1000000000, file=paste0(file, '_categories'),
			scaling=scaling, width=width, height=height, save_figure=TRUE)
	}
}


for (i in 1:length(file_paths)) {
	meta_list <- fread(
		paste0("~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/",
			file_root[i], "_figure_4.tsv.gz"))
	for (phe in unique(meta_list$phenotype)) {
		cat(phe, "\n")
		meta_list_tmp <- meta_list %>% filter(phenotype == phe)
		p <- make_manhattan_plot(meta_list_tmp$chromosome_name,
			meta_list_tmp$start_position,
			meta_list_tmp$Pvalue,
			threshold=1000, significance_T = 6.7e-7,
			label=meta_list_tmp$external_gene_name, 
			colour_1 = "#6583E6",
			colour_2 = "#384980")
		threshold <- ifelse(meta_list_tmp$case_control[1], 10, 10)
		p$p <- p$p + geom_label_repel(
			data = unique(subset(p$dt, y > threshold) %>% group_by(labels) %>% 
				filter(y == max(y))) %>% ungroup(),
			size = 2, aes(label=labels),
			color='grey30', box.padding = 0.2, force = 0.3,
			label.padding = 0.1, point.padding = 0.1, segment.color = 'grey50')
		width <- 230
		height <- 100
		scaling <- 1
		file <- paste0(phe, "_", file_root[i])
		ggsave(paste0(file, '.jpg'), p$p, width=width*scaling,
			height=height*scaling, units='mm')
		width <- 150
		print(meta_list_tmp %>% filter(Pvalue < 6.7e-7))
	}
}

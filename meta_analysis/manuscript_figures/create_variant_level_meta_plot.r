# Next, let's do the same thing, but for the variant level test statistics - Get to the end of this by today - this is my goal

files <- dir("~/Repositories/BRaVa_curation/meta_analysis/", full.names=TRUE)
variant_files <- grep(files, pattern="rare.exome.tsv.gz", value=TRUE)
variant_files <- setdiff(variant_files, "/Users/dpalmer/Repositories/BRaVa_curation/meta_analysis//MatHem_variant_meta_analysis1.rare.exome.tsv.gz")
dt_variant_results <- list()

for (file in variant_files) {
	cat(paste0(file, "\n"))
	dt_variant_results[[file]] <- fread(file, nrows = 10000) %>% mutate(file = file)
}
dt_variant_results <- rbindlist(dt_variant_results)
dt_variant_results <- dt_variant_results %>% mutate(
	chr = gsub("^chr([0-9X]+):.*", "\\1", MarkerName),
	pos = gsub("^chr[0-9X]+:([0-9]+):.*", "\\1", MarkerName)
)

p <- make_manhattan_plot(dt_variant_results$chr,
			dt_variant_results$pos,
			dt_variant_results$`P-value`,
			threshold=1000, significance_T = 8e-9,
			save_figure=TRUE,
			colour_1 = "#6583E6",
			colour_2 = "#384980"
			)

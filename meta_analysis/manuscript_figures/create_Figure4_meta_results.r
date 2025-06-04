# Create the manhattan plots
# What do we need for this?
# The region, the gene name, the minimum P-value among the hits

# Install and load biomaRt
# install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
library(ggplot2)
devtools::install_github("mkanai/rgsutil")
library(rgsutil)
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")
source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")

# Meta-analysis files for plotting should have been copied from BMRC to gcloud
# gsutil cp -r gene gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/
# gsutil mv 'gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/*gz gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/meta-analysis/

get_category_colors <- function(category=NULL)
{
	color_infectious <- "#822B20"
	color_neoplasms <- "#D26529"
	color_immune <- "#C6B844"
	color_metabolic <- "#999999"
	color_behavioral <- "#B17D89"
	color_nervous <- "#77216F"
	color_eye <- "#97C1A9"
	color_ear <- "#1F2D58"
	color_circulatory <- "#4F572E"
	color_respiratory <- "#B87726"
	color_digestive <- "#83BAE4"
	color_dermatologic <- "#1657C4"
	color_musculoskeletal <- "#5891E0"
	color_urinary <- "#F2C4C2"
	color_pregnancy <- "#F3C44E"
	color_other_factors <- "#4D7F36"

	category_colors <- list(
		Cardiovascular = color_circulatory,
		`Sense organs`= color_eye,
		Respiratory = color_respiratory,
		Neoplasms = color_neoplasms,
		Genitourinary = color_urinary,
		Muscloskeletal = color_musculoskeletal,
		Gastrointestinal = color_digestive,
		Dermatological = color_dermatologic,
		`Endocrine/Metabolic` = color_metabolic
		)

	if (!is.null(category)) {
		return(category_colors[category])
	}
	return(category_colors)
}

make_gene_manhattan_category_plot <- function(dt, buffer=100000000,
	chr_lengths=chr_lengths_38, significance_T=6.7e-7, ggplot_theme=theme_bw,
	save_figure=FALSE, file='file_out', scaling=1, width=230, height=100, 
	print_p=FALSE)
{  	
    contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
    start_end <- get_start_and_end(chr_lengths)
    dt_contigs <- data.frame(contig=contigs_, start=start_end$start,
    	end=start_end$end)

    # Now, we also want to shift by phenotype category
    if (!(all(c("Pvalue", "phenotype_category", "chromosome_name") %in% names(dt)))) {
    	cat("Pvalue, phenotype_category, or chromosome_name not present in the passed data.table\n")
    }

    if (!("position" %in% names(dt))) {
    	dt <- dt %>% mutate(position = (
    		start_position + (end_position - start_position)/2))
    }

    # Take the vector of phenotype categories
    # Order them alphabetically
    phenotype_categories <- sort(unique(dt$phenotype_category))
    print(phenotype_categories)
    phenotype_category_list <- list()
    i <- 1
   	for (phenotype_category in phenotype_categories) {
   		phenotype_category_list[[phenotype_category]] <- i
   		i <- i+1
   	}
    dt_categories <- data.frame(
    	phenotype_category = phenotype_categories,
    	start = ((dt_contigs$end[nrow(dt_contigs)]) * (seq(1,length(phenotype_categories))-1) + 1),
    	end = ((dt_contigs$end[nrow(dt_contigs)]) * (seq(1,length(phenotype_categories))))) %>%
    	mutate(middle = floor(start + (end-start)/2), length = (end-start)) %>% 
    	mutate(shifted_position = middle + (unlist(phenotype_category_list[phenotype_category])-1) * buffer)

    dt_plot <- dt %>% transmute(
    	contig=gsub("chr", "", chromosome_name),
    	position=as.integer(position),
    	pval=as.numeric(Pvalue),
    	labels=external_gene_name,
    	phenotype_category=phenotype_category) %>% mutate(
    	x = dt_categories[unlist(phenotype_category_list[phenotype_category]), 'start'] + 
    		dt_contigs[gsub('X', 23, contig), 'start'] + 
    		position + (unlist(phenotype_category_list[phenotype_category])-1) * buffer)

    p <- ggplot(dt_plot, aes(x=x, y=ifelse(pval < 1e-300, -300, -log10(pval)), col=phenotype_category)) + geom_point_rast(size=0.5)
    p <- p + geom_hline(yintercept=-log10(significance_T), color='#E15759', linetype='dashed') +
        scale_x_continuous(breaks=dt_categories$shifted_position, labels=dt_categories$phenotype_category) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        scale_color_manual(values = get_category_colors()) +
        labs(x='', y=TeX("$-\\log_{10}(P)$")) + ggplot_theme() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none")

    if (save_figure) {
        ggsave(paste0(file, '.jpg'), p, width=width*scaling, height=height*scaling, units='mm')
    }

    if (print_p) { print(p) }

    return(list(p=p, dt=dt_plot))
}

# Connect to Ensembl BioMart
ensembl <- useMart("ENSEMBL_MART_ENSEMBL")

# Choose the dataset for human genes
ensembl_dataset <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

# Retrieve Ensembl IDs, start and end positions
gene_info <- getBM(attributes = c(
	"ensembl_gene_id",
	"external_gene_name",
	"chromosome_name",
	"start_position",
	"end_position"),
	mart = ensembl_dataset)

# Display the first few rows of the result
gene_info <- data.table(gene_info, key = "ensembl_gene_id")

# Download the results files using Masa's rgsutil
file_paths <- c(
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/meta-analysis/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/AFR/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/AMR/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/EAS/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/EUR/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/SAS/",
	"gs://brava-meta-pilot-analysis/pilot-traits-may-2025-freeze-meta-analysis/gene/n_cases_100/non_EUR/"
)

file_root <- c("meta_analysis", "AFR", "AMR", "EAS", "EUR", "SAS", "non_EUR")

for (i in 1:length(file_paths)) {
	# First, grab the files
	files <- system(paste("gsutil ls", file_paths[i]), intern=TRUE)
	meta_list <- list()

	for (file in files) {
		cat(file, "\n")
		phenotype <- gsub(".*/(.*)_gene_meta_analysis_.*", "\\1", file)
		meta_list[[file]] <- rgsutil::read_gsfile(file) %>% filter(
			max_MAF %in% c("1e-04", "0.001"),
			Group %in% c(
				"damaging_missense_or_protein_altering",
				"pLoF",
				"pLoF;damaging_missense_or_protein_altering")) %>%
		filter(
			((class == "Burden") & (type == "Inverse variance weighted")) | 
			((class != "Burden") & (type == "Stouffer"))) %>%
		mutate(phenotype = phenotype) %>% group_by(Region, phenotype) %>% filter(Pvalue == min(Pvalue))
	}

	meta_list <- rbindlist(meta_list) %>% rename(ensembl_gene_id = Region)
	setkey(meta_list, "ensembl_gene_id")
	meta_list <- merge(gene_info, meta_list)
	meta_list <- meta_list %>% filter(chromosome_name %in% c(seq(1,22), "X"))

	# Remove the mosaic genes?
	meta_list <- meta_list %>% filter(!ensembl_gene_id %in% c("ENSG00000168769", "ENSG00000119772", "ENSG00000171456"))
	meta_list <- meta_list %>% mutate(Pvalue = ifelse(Pvalue == 0, 1e-320, Pvalue))
	meta_list <- meta_list %>% mutate(phenotype = gsub("_.*", "", phenotype))
	meta_list <- meta_list %>% mutate(
		case_control = ifelse(phenotype %in% case_ctrl, TRUE,
		ifelse(phenotype %in% cts, FALSE, NA)))

	# Write the information to disk, to speed up recreation of the plots.
	fwrite(meta_list, quote=FALSE, sep='\t',
		file=paste0("~/Repositories/BRaVa_curation/data/meta_analysis/meta_results/", file_root[i], "_figure_4.tsv.gz"))
}

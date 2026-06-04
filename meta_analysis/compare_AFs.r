library(data.table)
library(dplyr)
library(forcats)

source("~/Repositories/BRaVa_curation/QC/utils/pretty_plotting.r")
source("~/Repositories/BRaVa_curation/meta_analysis/meta_analysis_utils.r")

analysis_results_folder <- "~/Repositories/BRaVa_curation/data/meta_analysis/gene"
analysis_results_path_regexp <- "*cleaned*"
files <- dir(analysis_results_folder,
    pattern=analysis_results_path_regexp,
    full.names=TRUE)

dt_list <- list()
i <- 1
for (file in files) {
    file_info <- extract_file_info(gsub(".*/(.*)", "\\1", file))
    print(file_info)
    pheno <- file_info$phenotype
    dataset <- file_info$dataset
    ancestry <- file_info$ancestry
    n <- as.integer(file_info$n_controls) + as.integer(file_info$n_cases)
    dt_list[[i]] <- fread(cmd = paste("gzcat", file)) %>% 
    	select(Region, Group, max_MAF, MAC) %>% 
    	mutate(biobank = dataset, ancestry=ancestry, N=n) %>%
    	filter(Group != "Cauchy")
    i <- i+1
}

dt <- rbindlist(dt_list)
# Remove GEL, because they are unable to provide allele counts
dt <- dt %>% filter(biobank != "gel")
dt <- dt %>% filter(!grepl(";", Group))
dt <- dt %>% filter(!((biobank == "all-of-us") & !(Group %in% c("pLoF", "synonymous"))))
# In all of us, restrict to the synonymous and pLoFs for comparison
dt <- dt %>% mutate(CAF = (MAC/2)/N)

# Refactor the column
dt$Group <- recode(dt$Group, !!!renaming_plot_group_list)

# Print the updated data frame
pdf(file="compare_CAFs.pdf", width=10, height=10)
for (anc in unique(dt$ancestry)) {
	dt_tmp <- dt %>% filter(ancestry == anc)
	biobanks <- unique(dt_tmp$biobank)
	if (length(biobanks) > 1) {
		for (i in 1:(length(biobanks)-1)) {
			b1 <- biobanks[i]
			dt_1 <- dt_tmp %>% filter(
				biobank == b1,
				ancestry == anc) %>% 
				select(Region, max_MAF, CAF, Group)
			dt_1 <- data.table(dt_1)
			setkeyv(dt_1, c("Region", "max_MAF", "Group"))
			for (j in (i+1):length(biobanks)) {
				b2 <- biobanks[j]
				dt_2 <- dt_tmp %>% filter(
				biobank == b2,
				ancestry == anc) %>% 
				select(Region, max_MAF, CAF, Group)
				dt_2 <- data.table(dt_2)
				setkeyv(dt_2, c("Region", "max_MAF", "Group"))
				dt_plot <- merge(dt_1, dt_2)
				p <- ggplot(dt_plot, aes(x = `CAF.x`, y = `CAF.y`)) + 
				stat_bin_2d() + theme_classic() + 
				geom_abline(intercept=0, slope=1, col='indianred3') +
				labs(x=b1, y=b2, title=anc) + scale_fill_gradient(trans='log10', name='Frequency') +
				facet_wrap(vars(max_MAF, Group), scales = "free")
				print(p)
				p <- ggplot(dt_plot, aes(x = `CAF.x`, y = `CAF.y`)) + 
				stat_bin_2d() + theme_classic() + 
				geom_abline(intercept=0, slope=1, col='indianred3') +
				labs(x=b1, y=b2, title=anc) + scale_fill_gradient(trans='log10', name='Frequency') +
				facet_wrap(vars(max_MAF, Group), scales = "free") +
				scale_x_continuous(trans='log10') + 
				scale_y_continuous(trans='log10')
				print(p)
			}
		}
	}
}
dev.off()


# Print the updated data frame
dt <- dt %>% filter(max_MAF == 0.001, Group %in% c("pLoF", "Synonymous"))
pdf(file="compare_CAFs_0_001.pdf", width=10, height=5)
for (anc in unique(dt$ancestry)) {
	dt_tmp <- dt %>% filter(ancestry == anc)
	biobanks <- unique(dt_tmp$biobank)
	if (length(biobanks) > 1) {
		for (i in 1:(length(biobanks)-1)) {
			b1 <- biobanks[i]
			dt_1 <- dt_tmp %>% filter(
				biobank == b1,
				ancestry == anc) %>% 
				select(Region, max_MAF, CAF, Group)
			dt_1 <- data.table(dt_1)
			setkeyv(dt_1, c("Region", "max_MAF", "Group"))
			for (j in (i+1):length(biobanks)) {
				b2 <- biobanks[j]
				dt_2 <- dt_tmp %>% filter(
				biobank == b2,
				ancestry == anc) %>% 
				select(Region, max_MAF, CAF, Group)
				dt_2 <- data.table(dt_2)
				setkeyv(dt_2, c("Region", "max_MAF", "Group"))
				dt_plot <- merge(dt_1, dt_2)
				p <- ggplot(dt_plot, aes(x = `CAF.x`, y = `CAF.y`)) + 
				stat_bin_2d() + theme_classic() + 
				geom_abline(intercept=0, slope=1, col='indianred3') +
				labs(x=b1, y=b2, title=anc) + scale_fill_gradient(trans='log10', name='Frequency') +
				facet_wrap(vars(max_MAF, Group), scales = "free")
				print(p)
				p <- ggplot(dt_plot, aes(x = `CAF.x`, y = `CAF.y`)) + 
				stat_bin_2d() + theme_classic() + 
				geom_abline(intercept=0, slope=1, col='indianred3') +
				labs(x=b1, y=b2, title=anc) + scale_fill_gradient(trans='log10', name='Frequency') +
				facet_wrap(vars(max_MAF, Group), scales = "free") +
				scale_x_continuous(trans='log10') + 
				scale_y_continuous(trans='log10')
				print(p)
			}
		}
	}
}
dev.off()


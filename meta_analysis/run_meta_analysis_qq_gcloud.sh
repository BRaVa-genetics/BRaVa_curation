#!/bin/bash

# Rscript meta_analysis_qq.r --meta_analysis_results_folder ../data/meta_analysis/gcloud/meta_results/n_cases_100 --out gcloud_meta_analysis_qq_plots_100.pdf
# Rscript meta_analysis_qq.r --meta_analysis_results_folder ../data/meta_analysis/gcloud/meta_results/n_cases_1000 --out gcloud_meta_analysis_qq_plots_1000.pdf

Rscript meta_analysis_qq.r --type 'weighted Fisher' --meta_analysis_results_folder ../data/meta_analysis/meta_results/n_cases_100 --out gcloud_meta_analysis_qq_plots_100_new.pdf
# Rscript meta_analysis_qq.r --type 'weighted Fisher' --meta_analysis_results_folder ../data/meta_analysis/gcloud/meta_results/n_cases_1000 --out gcloud_meta_analysis_qq_plots_1000.pdf

# Plot each biobank's results
Rscript analysis_qq.r --analysis_results_folder ../data/meta_analysis/gcloud --analysis_results_path_regexp "*cleaned*" --out gcloud_biobanks_qq_plots.pdf

#!/usr/bin/env bash
#SBATCH --job-name=meta-analysis-qq
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

META_FILE_GENE=$1
OUT=$2

echo "Rscript meta_analysis_qq.r --meta_analysis_results_file $META_FILE_GENE --include_gene_names --out $OUT"
Rscript meta_analysis_qq.r --meta_analysis_results_file $META_FILE_GENE --out $OUT --include_gene_names

new_OUT=$(echo "$OUT" | sed "s/\.pdf/_burden.pdf/")

echo "Rscript meta_analysis_qq.r --meta_analysis_results_file $META_FILE_GENE --include_gene_names --burden_only_plot --out $new_OUT"
Rscript meta_analysis_qq.r --meta_analysis_results_file $META_FILE_GENE --out $new_OUT --include_gene_names --burden_only_plot

# new_META_FILE_GENE=$(echo $META_FILE_GENE | sed "s/_cutoff.tsv.gz/_cutoff_extra_cauchy.tsv.gz/g")
# new_OUT=$(echo "$OUT" | sed "s/\.pdf/_cauchy.pdf/")
# echo "Rscript meta_analysis_cauchy_qq.r --meta_analysis_results_file $META_FILE_GENE --include_gene_names --out $new_OUT"
# Rscript meta_analysis_cauchy_qq.r --meta_analysis_results_file $new_META_FILE_GENE --out $new_OUT --include_gene_names

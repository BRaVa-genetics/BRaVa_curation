#!/usr/bin/env bash
#SBATCH --job-name=analysis-qq
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

FILE_GENE=$1
OUT=$2
new_OUT=$(echo "$OUT" | sed "s/\.pdf/_burden.pdf/")

echo "Rscript analysis_qq.r --analysis_results_file $FILE_GENE --include_gene_names --out $OUT"
Rscript analysis_qq.r --analysis_results_file $FILE_GENE --include_gene_names --out $OUT

echo "Rscript analysis_qq.r --analysis_results_file $FILE_GENE --include_gene_names --out $new_OUT"
Rscript analysis_qq.r --analysis_results_file $FILE_GENE --burden_only_plot --include_gene_names --out $new_OUT

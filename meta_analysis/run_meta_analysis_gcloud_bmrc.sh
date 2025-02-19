#!/usr/bin/env bash
#SBATCH --job-name=meta-analysis
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

module load R

FILES_GENE=$1
OUT=$2

Rscript meta_analysis.r --file_paths $FILES_GENE --no_sex_check --out $OUT

# new_FILES_GENE=$(echo "$FILES_GENE" | sed "s/txt.gz/extra_cauchy.gz/g")
# new_OUT=$(echo $OUT | sed "s/cutoff/cutoff_extra_cauchy/g")

# Rscript meta_analysis.r --file_paths $new_FILES_GENE --no_sex_check --out $new_OUT

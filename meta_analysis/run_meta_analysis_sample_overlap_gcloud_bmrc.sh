#!/usr/bin/env bash
#SBATCH --job-name=meta-analysis
#SBATCH --output=/well/lindgren/dpalmer/logs/%x_%j.log
#SBATCH --error=/well/lindgren/dpalmer/logs/%x_%j.err

module load R

FILES_GENE=$1
OUT=$2

Rscript meta_analysis_sample_overlap.r --file_paths $FILES_GENE --no_sex_check --out $OUT
